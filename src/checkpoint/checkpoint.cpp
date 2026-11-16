#include "checkpoint/checkpoint.hpp"

#include <atomic>
#include <bit>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>

#if !defined(_WIN32)
#include <csignal>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>
#else
#include <io.h>
#endif

namespace qraft::internal {
namespace {

constexpr const char* kMagic = "QRAFT-CHECKPOINT 1";

/// A journaled bipartition: the positions, among the subproblem's ascending
/// taxa, of its smaller side, and whether that side is the left one.
struct SplitEntry {
  int active = 0;
  bool listed_left = true;
  std::vector<int> positions;
};

/// Process-wide state of one journaling run. Runs that journal are exclusive,
/// so one run owns it at a time. The replay maps are filled before the run
/// starts and only read afterwards.
struct CheckpointState {
  std::mutex mutex;
  std::condition_variable wake;
  bool active = false;
  bool journaling = false;
  bool finished = false;
  std::string path;
  std::string progress_file;
  double interval_s = 1800.0;
  double progress_interval_s = 30.0;
  std::string header;

  // Decisions made by this run and not yet written, and every decision the
  // file holds (the resumed run's journal first).
  std::vector<std::string> pending;
  std::vector<std::string> written;
  bool flush_requested = false;
  bool stop = false;
  std::thread monitor;

  std::chrono::steady_clock::time_point started;
  std::chrono::steady_clock::time_point last_save;
  std::atomic<long long> last_save_unix{0};
  std::atomic<long long> saves{0};
  std::atomic<long long> recorded{0};
  std::atomic<long long> saved_decisions{0};
  std::atomic<long long> replayed{0};
  long long loaded_decisions = 0;
  bool resumed = false;
  std::atomic<int> stage{0};
  std::atomic<int> invocations[3] = {0, 0, 0};
  bool write_failed = false;
  // The advisory lock on path.lock that keeps two runs off one journal.
  int lock_fd = -1;

  // Test hooks: QRAFT_CHECKPOINT_TEST_KILL_AFTER=N kills the process after its
  // N-th journal write, QRAFT_CHECKPOINT_TEST_KILL_MID=N during it (after the
  // temporary file is complete, before the rename).
  long long kill_after = 0;
  long long kill_mid = 0;

  std::map<std::pair<int, uint64_t>, SplitEntry> splits;
  std::map<std::pair<int, uint64_t>, std::array<int, 4>> regrafts;
  std::map<int, std::vector<std::vector<NniMove>>> nni;
  std::map<std::pair<int, int>, SprRoundLog> spr;
};

CheckpointState& state() {
  static CheckpointState instance;
  return instance;
}

const char* stage_name(int stage) {
  switch (static_cast<RunStage>(stage)) {
    case RunStage::Startup: return "startup";
    case RunStage::StartingTree: return "starting-tree";
    case RunStage::Divide: return "divide-regraft";
    case RunStage::Nni: return "nni";
    case RunStage::Spr: return "spr";
    case RunStage::Support: return "support";
    case RunStage::Done: return "done";
  }
  return "unknown";
}

uint64_t fnv1a(uint64_t hash, const std::string& text) {
  for (const unsigned char c : text) {
    hash ^= c;
    hash *= 1099511628211ULL;
  }
  hash ^= 0x0a;
  hash *= 1099511628211ULL;
  return hash;
}

constexpr uint64_t kFnvOffset = 14695981039346656037ULL;

long long env_count(const char* name) {
  const char* value = std::getenv(name);
  if (value == nullptr || *value == '\0') return 0;
  return std::strtoll(value, nullptr, 10);
}

[[noreturn]] void kill_self() {
#if !defined(_WIN32)
  std::raise(SIGKILL);
#endif
  std::_Exit(137);
}

std::string json_escape(const std::string& text) {
  std::string out;
  out.reserve(text.size());
  for (const char c : text) {
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          char buffer[8];
          std::snprintf(buffer, sizeof(buffer), "\\u%04x",
                        static_cast<unsigned>(static_cast<unsigned char>(c)));
          out += buffer;
        } else {
          out += c;
        }
    }
  }
  return out;
}

/// The algorithmic options a resumed run must share with its journal. The
/// thread count, progress and output paths may differ.
std::string options_fingerprint(const Options& options) {
  std::ostringstream out;
  out << std::setprecision(17);
  out << "divide=" << options.divide
      << " divide_score=" << static_cast<int>(options.divide_score)
      << " test=" << options.bipartition_test.enabled << ','
      << options.bipartition_test.coefficient << ','
      << options.bipartition_test.crossover << ','
      << options.bipartition_test.proposals << ','
      << options.bipartition_test.anchor_pairs
      << " regraft_score=" << static_cast<int>(options.regraft_score)
      << " nni=" << options.nni << ',' << options.nni_move_cap << ','
      << options.nni_support_gate << ',' << options.nni_batch
      << " spr=" << options.spr_refine << ',' << options.spr_reads << ','
      << options.spr_min_clade << ',' << options.spr_candidates << ','
      << options.spr_rounds
      << " refine_start=" << options.refine_start
      << " support=" << static_cast<int>(options.support) << ','
      << static_cast<int>(options.support_frequencies) << ','
      << options.support_detail
      << " starting_tree=" << options.starting_tree_path
      << " guide_tree=" << options.guide_tree_path;
  return out.str();
}

std::string make_header(const Options& options,
                        const std::string& input_fingerprint) {
  std::ostringstream out;
  out << kMagic << '\n';
  out << "version " << options.program_version << '\n';
  out << "options " << options_fingerprint(options) << '\n';
  out << "input " << input_fingerprint << '\n';
  return out.str();
}

#if !defined(_WIN32)
void fsync_directory_of(const std::string& path) {
  std::filesystem::path dir = std::filesystem::path(path).parent_path();
  if (dir.empty()) dir = ".";
  const int fd = ::open(dir.c_str(), O_RDONLY);
  if (fd < 0) return;
  ::fsync(fd);
  ::close(fd);
}
#endif

/// Writes the whole journal to path.tmp, flushes it to disk and renames it
/// over `path`, so a cut leaves either the previous journal or this one.
void write_journal_file(CheckpointState& s,
                        const std::vector<std::string>& lines) {
  const std::string tmp = s.path + ".tmp";
  std::FILE* file = std::fopen(tmp.c_str(), "wb");
  if (file == nullptr) {
    throw std::runtime_error("cannot write checkpoint journal: " + tmp);
  }
  uint64_t hash = kFnvOffset;
  bool ok = std::fwrite(s.header.data(), 1, s.header.size(), file) ==
            s.header.size();
  for (const std::string& line : lines) {
    if (!ok) break;
    hash = fnv1a(hash, line);
    ok = std::fwrite(line.data(), 1, line.size(), file) == line.size() &&
         std::fputc('\n', file) != EOF;
  }
  if (ok) {
    ok = std::fprintf(file, "END %zu %016llx\n", lines.size(),
                      static_cast<unsigned long long>(hash)) > 0;
  }
  ok = ok && std::fflush(file) == 0;
#if !defined(_WIN32)
  ok = ok && ::fsync(::fileno(file)) == 0;
#else
  ok = ok && ::_commit(::_fileno(file)) == 0;
#endif
  ok = (std::fclose(file) == 0) && ok;
  if (!ok) {
    throw std::runtime_error("cannot write checkpoint journal: " + tmp);
  }
  const long long save_number = s.saves.load(std::memory_order_relaxed) + 1;
  if (s.kill_mid > 0 && save_number == s.kill_mid) kill_self();
  std::error_code error;
  std::filesystem::rename(tmp, s.path, error);
  if (error) {
    throw std::runtime_error("cannot replace checkpoint journal " + s.path +
                             ": " + error.message());
  }
#if !defined(_WIN32)
  fsync_directory_of(s.path);
#endif
  s.saves.store(save_number, std::memory_order_relaxed);
  s.saved_decisions.store(static_cast<long long>(lines.size()),
                          std::memory_order_relaxed);
  s.last_save_unix.store(
      std::chrono::duration_cast<std::chrono::seconds>(
          std::chrono::system_clock::now().time_since_epoch())
          .count(),
      std::memory_order_relaxed);
  if (s.kill_after > 0 && save_number == s.kill_after) kill_self();
}

/// Moves the pending decisions into the journal and writes it. Called with
/// the lock held; writes without it. Only the monitor thread, or the owner
/// after the monitor has stopped, calls it.
void flush_locked(CheckpointState& s, std::unique_lock<std::mutex>& lock) {
  for (std::string& line : s.pending) s.written.push_back(std::move(line));
  s.pending.clear();
  lock.unlock();
  try {
    write_journal_file(s, s.written);
  } catch (const std::exception& error) {
    if (!s.write_failed) {
      std::cerr << "qraft: warning: " << error.what() << '\n';
    }
    s.write_failed = true;
  }
  lock.lock();
  s.last_save = std::chrono::steady_clock::now();
  s.flush_requested = false;
}

long long read_rss_bytes(long long* peak) {
  *peak = 0;
#if defined(__linux__)
  long long resident = 0;
  {
    std::ifstream statm("/proc/self/statm");
    long long size_pages = 0, resident_pages = 0;
    if (statm >> size_pages >> resident_pages) {
      resident = resident_pages * static_cast<long long>(::sysconf(_SC_PAGESIZE));
    }
  }
  std::ifstream status("/proc/self/status");
  std::string line;
  while (std::getline(status, line)) {
    if (line.rfind("VmHWM:", 0) == 0) {
      *peak = std::strtoll(line.c_str() + 6, nullptr, 10) * 1024;
      break;
    }
  }
  return resident;
#else
  return 0;
#endif
}

/// Writes the status file by atomic replace; a reader never sees half of it.
void write_status(CheckpointState& s, const char* run_state) {
  if (s.progress_file.empty()) return;
  RunStats& stats = run_stats();
  const auto load = [](const std::atomic<long long>& counter) {
    return counter.load(std::memory_order_relaxed);
  };
  long long peak = 0;
  const long long rss = read_rss_bytes(&peak);
  const double elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    s.started)
          .count();
  std::ostringstream out;
  out << "{\n";
  out << "  \"state\": \"" << run_state << "\",\n";
  out << "  \"stage\": \"" << stage_name(s.stage.load(std::memory_order_relaxed))
      << "\",\n";
  out << "  \"elapsed_s\": " << std::fixed << std::setprecision(1) << elapsed
      << ",\n";
  out << "  \"subproblems\": " << load(stats.subproblems) << ",\n";
  out << "  \"completed_subproblems\": " << load(stats.completed_subproblems)
      << ",\n";
  out << "  \"regraft_calls\": " << load(stats.regraft_calls) << ",\n";
  out << "  \"max_depth\": " << load(stats.max_depth) << ",\n";
  out << "  \"nni_moves\": " << load(stats.nni_moves_accepted) << ",\n";
  out << "  \"spr_rounds\": " << load(stats.spr_rounds) << ",\n";
  out << "  \"spr_candidates\": " << load(stats.spr_candidates_scored) << ",\n";
  out << "  \"spr_moves\": " << load(stats.spr_moves) << ",\n";
  out << "  \"rss_bytes\": " << rss << ",\n";
  out << "  \"peak_rss_bytes\": " << peak << ",\n";
  out << "  \"checkpoint\": {";
  if (s.journaling) {
    out << "\"path\": \"" << json_escape(s.path) << "\", "
        << "\"resumed\": " << (s.resumed ? "true" : "false") << ", "
        << "\"loaded_decisions\": " << s.loaded_decisions << ", "
        << "\"replayed\": " << s.replayed.load(std::memory_order_relaxed) << ", "
        << "\"decisions\": "
        << s.loaded_decisions + s.recorded.load(std::memory_order_relaxed)
        << ", "
        << "\"saved_decisions\": "
        << s.saved_decisions.load(std::memory_order_relaxed) << ", "
        << "\"saves\": " << s.saves.load(std::memory_order_relaxed) << ", "
        << "\"last_save_unix\": "
        << s.last_save_unix.load(std::memory_order_relaxed) << ", "
        << "\"interval_s\": " << s.interval_s;
  }
  out << "},\n";
  out << "  \"updated_unix\": "
      << std::chrono::duration_cast<std::chrono::seconds>(
             std::chrono::system_clock::now().time_since_epoch())
             .count()
      << "\n}\n";
  const std::string tmp = s.progress_file + ".tmp";
  {
    std::ofstream file(tmp, std::ios::binary | std::ios::trunc);
    if (!file) return;
    file << out.str();
    if (!file) return;
  }
  std::error_code error;
  std::filesystem::rename(tmp, s.progress_file, error);
}

void monitor_loop(CheckpointState& s) {
  std::unique_lock<std::mutex> lock(s.mutex);
  auto next_status = std::chrono::steady_clock::now();
  const auto status_period = std::chrono::duration_cast<
      std::chrono::steady_clock::duration>(
      std::chrono::duration<double>(s.progress_interval_s));
  while (!s.stop) {
    if (s.flush_requested) {
      flush_locked(s, lock);
      continue;
    }
    const auto now = std::chrono::steady_clock::now();
    if (!s.progress_file.empty() && now >= next_status) {
      lock.unlock();
      write_status(s, "running");
      lock.lock();
      next_status = now + status_period;
      continue;
    }
    if (s.progress_file.empty()) {
      s.wake.wait(lock, [&] { return s.stop || s.flush_requested; });
    } else {
      s.wake.wait_until(lock, next_status,
                        [&] { return s.stop || s.flush_requested; });
    }
  }
}

void record_line(std::string line) {
  CheckpointState& s = state();
  std::lock_guard<std::mutex> lock(s.mutex);
  s.pending.push_back(std::move(line));
  s.recorded.fetch_add(1, std::memory_order_relaxed);
  if (!s.flush_requested) {
    const double since =
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      s.last_save)
            .count();
    if (since >= s.interval_s) {
      s.flush_requested = true;
      s.wake.notify_one();
    }
  }
}

/// Takes an exclusive advisory lock on path.lock for the run, so that two
/// runs never write one journal. Windows relies on the user (no lock).
void lock_journal(CheckpointState& s, const std::string& path) {
#if !defined(_WIN32)
  const std::string lock_path = path + ".lock";
  const int fd = ::open(lock_path.c_str(), O_RDWR | O_CREAT, 0644);
  if (fd < 0) {
    throw std::runtime_error("cannot create checkpoint lock file: " + lock_path);
  }
  if (::flock(fd, LOCK_EX | LOCK_NB) != 0) {
    ::close(fd);
    throw std::runtime_error("another qraft run is using checkpoint journal " +
                             path);
  }
  s.lock_fd = fd;
#else
  (void)s;
  (void)path;
#endif
}

void unlock_journal(CheckpointState& s) {
#if !defined(_WIN32)
  if (s.lock_fd >= 0) ::close(s.lock_fd);
#endif
  s.lock_fd = -1;
}

void open_run(CheckpointState& s, const Options& options,
              const std::string& input_fingerprint);

[[noreturn]] void bad_journal(const std::string& path, const std::string& why) {
  throw std::runtime_error("checkpoint journal " + path + " is not usable: " +
                           why);
}

/// Reads a journal, checks it belongs to this run, and fills the replay maps.
void load_journal(CheckpointState& s) {
  std::ifstream in(s.path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("cannot read checkpoint journal: " + s.path);
  }
  std::string line;
  std::istringstream expected(s.header);
  std::string expected_line;
  while (std::getline(expected, expected_line)) {
    if (!std::getline(in, line)) bad_journal(s.path, "it ends in its header");
    if (line != expected_line) {
      const std::string field = expected_line.substr(0, expected_line.find(' '));
      bad_journal(s.path, "its " + field +
                              " line differs from this run's (journal: \"" +
                              line + "\", this run: \"" + expected_line +
                              "\")");
    }
  }
  uint64_t hash = kFnvOffset;
  std::vector<std::string> lines;
  bool ended = false;
  while (std::getline(in, line)) {
    if (line.rfind("END ", 0) == 0) {
      std::istringstream end(line.substr(4));
      std::size_t count = 0;
      std::string hex;
      end >> count >> hex;
      std::ostringstream want;
      want << std::hex << std::setw(16) << std::setfill('0') << hash;
      if (count != lines.size() || hex != want.str()) {
        bad_journal(s.path, "its checksum does not match");
      }
      ended = true;
      break;
    }
    hash = fnv1a(hash, line);
    lines.push_back(line);
  }
  if (!ended) bad_journal(s.path, "it has no END line");

  for (const std::string& record : lines) {
    std::istringstream fields(record);
    char kind = 0;
    fields >> kind;
    if (kind == 'S') {
      int inv = 0, active = 0, listed_left = 0, count = 0;
      uint64_t path = 0;
      fields >> inv >> path >> active >> listed_left >> count;
      SplitEntry entry;
      entry.active = active;
      entry.listed_left = listed_left != 0;
      entry.positions.resize(static_cast<std::size_t>(std::max(0, count)));
      for (int& p : entry.positions) fields >> p;
      if (!fields || !s.splits.emplace(std::make_pair(inv, path), entry).second) {
        bad_journal(s.path, "a bipartition record is malformed or repeated");
      }
    } else if (kind == 'G') {
      int inv = 0;
      uint64_t path = 0;
      std::array<int, 4> joins{};
      fields >> inv >> path >> joins[0] >> joins[1] >> joins[2] >> joins[3];
      if (!fields ||
          !s.regrafts.emplace(std::make_pair(inv, path), joins).second) {
        bad_journal(s.path, "a regraft record is malformed or repeated");
      }
    } else if (kind == 'N') {
      int inv = 0, count = 0;
      fields >> inv >> count;
      std::vector<NniMove> group(static_cast<std::size_t>(std::max(0, count)));
      for (NniMove& move : group) {
        uint64_t gain_bits = 0;
        fields >> move.edge_a >> move.edge_b >> move.swap_a_neighbor >>
            move.swap_b_neighbor >> move.alternative >> gain_bits;
        move.gain = std::bit_cast<double>(gain_bits);
      }
      if (!fields) bad_journal(s.path, "an NNI record is malformed");
      s.nni[inv].push_back(std::move(group));
    } else if (kind == 'P') {
      int inv = 0, round = 0, moved = 0;
      SprOutcome outcome;
      fields >> inv >> round >> outcome.index >> moved >> outcome.edge_a >>
          outcome.edge_b >> outcome.joins[0] >> outcome.joins[1] >>
          outcome.joins[2] >> outcome.joins[3];
      outcome.moved = moved != 0;
      SprRoundLog& log = s.spr[{inv, round}];
      if (!fields || log.complete ||
          (!log.outcomes.empty() && log.outcomes.back().index >= outcome.index)) {
        bad_journal(s.path, "an SPR record is malformed or out of order");
      }
      log.outcomes.push_back(outcome);
    } else if (kind == 'Q') {
      int inv = 0, round = 0;
      fields >> inv >> round;
      if (!fields) bad_journal(s.path, "an SPR round record is malformed");
      s.spr[{inv, round}].complete = true;
    } else if (kind == 'D') {
      // The journaled run finished; its decisions replay to the same tree.
    } else {
      bad_journal(s.path, "it holds an unknown record");
    }
  }
  s.loaded_decisions = static_cast<long long>(lines.size());
  s.written = std::move(lines);
}

}  // namespace

void checkpoint_note_stage(RunStage stage) {
  state().stage.store(static_cast<int>(stage), std::memory_order_relaxed);
}

uint64_t checkpoint_genes_hash(const std::vector<GeneTree>& genes) {
  uint64_t hash = kFnvOffset;
  for (const GeneTree& gene : genes) {
    hash = fnv1a(hash, "tree " + std::to_string(gene.root));
    for (const GeneTree::Node& node : gene.nodes) {
      std::string text = node.label;
      for (int neighbor : node.neighbors) {
        text += ' ';
        text += std::to_string(neighbor);
      }
      hash = fnv1a(hash, text);
    }
  }
  return hash;
}

std::string checkpoint_input_fingerprint(const Options& options,
                                         const TaxonRegistry& taxa,
                                         std::size_t gene_count,
                                         std::size_t physical_nodes,
                                         uint64_t genes_hash) {
  if (options.checkpoint_path.empty()) return std::string();
  uint64_t labels = kFnvOffset;
  for (const std::string& label : taxa.labels()) labels = fnv1a(labels, label);
  // The starting and guide trees are identified by their content, not only
  // by their path.
  const auto file_hash = [](const std::string& path) {
    if (path.empty()) return std::string("-");
    std::ifstream in(path, std::ios::binary);
    if (!in) return std::string("unreadable");
    uint64_t hash = kFnvOffset;
    std::string line;
    while (std::getline(in, line)) hash = fnv1a(hash, line);
    std::ostringstream out;
    out << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
  };
  std::ostringstream out;
  out << "genes=" << gene_count << " taxa=" << taxa.size()
      << " nodes=" << physical_nodes << std::hex << std::setfill('0')
      << " gene_trees=" << std::setw(16) << genes_hash
      << " labels=" << std::setw(16) << labels
      << " starting_tree=" << file_hash(options.starting_tree_path)
      << " guide_tree=" << file_hash(options.guide_tree_path);
  return out.str();
}

ScopedCheckpointRun::ScopedCheckpointRun(const Options& options,
                                         std::string input_fingerprint) {
  if (options.checkpoint_path.empty() && options.progress_file.empty()) {
    if (options.resume) {
      throw std::invalid_argument("--resume needs a checkpoint journal");
    }
    return;
  }
  if (options.resume && options.checkpoint_path.empty()) {
    throw std::invalid_argument("--resume needs a checkpoint journal");
  }
  if (!std::isfinite(options.checkpoint_interval_s) ||
      options.checkpoint_interval_s < 0.0) {
    throw std::invalid_argument(
        "the checkpoint interval must be a finite number of seconds >= 0");
  }
  if (!std::isfinite(options.progress_interval_s) ||
      options.progress_interval_s <= 0.0) {
    throw std::invalid_argument(
        "the progress interval must be a finite number of seconds > 0");
  }
  if (!options.checkpoint_path.empty() && !options.progress_file.empty()) {
    // The status file is replaced every few seconds; it must not be the
    // journal or one of its companion files.
    const auto canonical = [](const std::string& path) {
      std::error_code error;
      std::filesystem::path resolved = std::filesystem::weakly_canonical(path, error);
      return error ? std::filesystem::path(path) : resolved;
    };
    const std::string& journal = options.checkpoint_path;
    const std::string& status = options.progress_file;
    for (const std::string& a : {journal, journal + ".tmp", journal + ".lock"}) {
      for (const std::string& b : {status, status + ".tmp"}) {
        if (canonical(a) == canonical(b)) {
          throw std::invalid_argument(
              "the progress file must differ from the checkpoint journal and "
              "its .tmp and .lock files");
        }
      }
    }
  }
  CheckpointState& s = state();
  std::lock_guard<std::mutex> lock(s.mutex);
  if (s.active) throw std::logic_error("a checkpointed run is already active");
  if (!options.checkpoint_path.empty()) {
    lock_journal(s, options.checkpoint_path);
  }
  try {
    open_run(s, options, input_fingerprint);
  } catch (...) {
    unlock_journal(s);
    throw;
  }
  s.monitor = std::thread([&s] { monitor_loop(s); });
  active_ = true;
}

namespace {

void open_run(CheckpointState& s, const Options& options,
              const std::string& input_fingerprint) {
  s.journaling = !options.checkpoint_path.empty();
  s.finished = false;
  s.path = options.checkpoint_path;
  s.progress_file = options.progress_file;
  s.interval_s = options.checkpoint_interval_s;
  s.progress_interval_s = options.progress_interval_s;
  s.header = make_header(options, input_fingerprint);
  s.pending.clear();
  s.written.clear();
  s.flush_requested = false;
  s.stop = false;
  s.saves = 0;
  s.recorded = 0;
  s.saved_decisions = 0;
  s.replayed = 0;
  s.last_save_unix = 0;
  s.loaded_decisions = 0;
  s.resumed = false;
  s.write_failed = false;
  for (auto& counter : s.invocations) counter = 0;
  s.splits.clear();
  s.regrafts.clear();
  s.nni.clear();
  s.spr.clear();
  s.kill_after = env_count("QRAFT_CHECKPOINT_TEST_KILL_AFTER");
  s.kill_mid = env_count("QRAFT_CHECKPOINT_TEST_KILL_MID");
  if (s.journaling) {
    std::error_code error;
    const bool exists = std::filesystem::exists(s.path, error);
    if (options.resume && exists) {
      load_journal(s);
      s.resumed = true;
      s.saved_decisions = s.loaded_decisions;
    } else if (options.resume) {
      // A run cut before its first write left no journal: nothing finished
      // was lost beyond one interval, so the run starts from the beginning.
      std::cerr << "qraft: warning: no checkpoint journal at " << s.path
                << "; starting from the beginning\n";
    } else if (exists) {
      throw std::runtime_error(
          "checkpoint journal already exists: " + s.path +
          " (pass --resume to continue that run, or remove the file)");
    }
  }
  s.started = std::chrono::steady_clock::now();
  s.last_save = s.started;
  s.active = true;
}

}  // namespace

void ScopedCheckpointRun::finish() {
  if (!active_) return;
  CheckpointState& s = state();
  if (s.journaling) record_line("D");
  checkpoint_note_stage(RunStage::Done);
  {
    std::lock_guard<std::mutex> lock(s.mutex);
    s.finished = true;
  }
}

ScopedCheckpointRun::~ScopedCheckpointRun() {
  if (!active_) return;
  CheckpointState& s = state();
  {
    std::lock_guard<std::mutex> lock(s.mutex);
    s.stop = true;
  }
  s.wake.notify_all();
  if (s.monitor.joinable()) s.monitor.join();
  std::unique_lock<std::mutex> lock(s.mutex);
  // Completed work is written on every exit, finished or failed.
  if (s.journaling && !s.pending.empty()) flush_locked(s, lock);
  const bool finished = s.finished;
  lock.unlock();
  try {
    write_status(s, finished ? "done" : "failed");
  } catch (...) {
  }
  lock.lock();
  s.active = false;
  s.journaling = false;
  s.splits.clear();
  s.regrafts.clear();
  s.nni.clear();
  s.spr.clear();
  s.written.clear();
  unlock_journal(s);
}

bool checkpoint_journaling() {
  CheckpointState& s = state();
  return s.active && s.journaling;
}

int checkpoint_next_invocation(char kind) {
  if (!checkpoint_journaling()) return -1;
  const int slot = kind == 'R' ? 0 : kind == 'N' ? 1 : 2;
  return state().invocations[slot].fetch_add(1, std::memory_order_relaxed);
}

std::optional<Bitset> checkpoint_split(int invocation, uint64_t path,
                                       int active) {
  if (invocation < 0) return std::nullopt;
  const CheckpointState& s = state();
  const auto it = s.splits.find({invocation, path});
  if (it == s.splits.end()) return std::nullopt;
  const SplitEntry& entry = it->second;
  if (entry.active != active) {
    throw std::logic_error(
        "checkpoint replay: a journaled bipartition does not match its "
        "subproblem");
  }
  Bitset side(active);
  if (entry.listed_left) {
    for (int p : entry.positions) side.set(p);
  } else {
    std::vector<char> right(static_cast<std::size_t>(active), 0);
    for (int p : entry.positions) {
      if (p < 0 || p >= active) {
        throw std::logic_error("checkpoint replay: a bipartition is out of range");
      }
      right[static_cast<std::size_t>(p)] = 1;
    }
    for (int i = 0; i < active; ++i) {
      if (!right[static_cast<std::size_t>(i)]) side.set(i);
    }
  }
  return side;
}

void checkpoint_record_split(int invocation, uint64_t path,
                             const Bitset& left_side) {
  if (invocation < 0) return;
  const int active = left_side.size();
  const int left_count = left_side.count();
  const bool listed_left = left_count <= active - left_count;
  std::ostringstream out;
  out << "S " << invocation << ' ' << path << ' ' << active << ' '
      << (listed_left ? 1 : 0) << ' '
      << (listed_left ? left_count : active - left_count);
  for (int i = 0; i < active; ++i) {
    if (left_side.test(i) == listed_left) out << ' ' << i;
  }
  record_line(out.str());
}

std::optional<std::pair<RegraftEdge, RegraftEdge>> checkpoint_regraft(
    int invocation, uint64_t path) {
  if (invocation < 0) return std::nullopt;
  const CheckpointState& s = state();
  const auto it = s.regrafts.find({invocation, path});
  if (it == s.regrafts.end()) return std::nullopt;
  RegraftEdge left, right;
  left.a = it->second[0];
  left.b = it->second[1];
  right.a = it->second[2];
  right.b = it->second[3];
  return std::make_pair(left, right);
}

void checkpoint_record_regraft(int invocation, uint64_t path,
                               const RegraftEdge& left,
                               const RegraftEdge& right) {
  if (invocation < 0) return;
  std::ostringstream out;
  out << "G " << invocation << ' ' << path << ' ' << left.a << ' ' << left.b
      << ' ' << right.a << ' ' << right.b;
  record_line(out.str());
}

const std::vector<std::vector<NniMove>>* checkpoint_nni_groups(int invocation) {
  if (invocation < 0) return nullptr;
  const CheckpointState& s = state();
  const auto it = s.nni.find(invocation);
  return it == s.nni.end() ? nullptr : &it->second;
}

void checkpoint_record_nni_group(int invocation,
                                 const std::vector<NniMove>& group) {
  if (invocation < 0 || group.empty()) return;
  std::ostringstream out;
  out << "N " << invocation << ' ' << group.size();
  for (const NniMove& move : group) {
    out << ' ' << move.edge_a << ' ' << move.edge_b << ' '
        << move.swap_a_neighbor << ' ' << move.swap_b_neighbor << ' '
        << move.alternative << ' ' << std::bit_cast<uint64_t>(move.gain);
  }
  record_line(out.str());
}

const SprRoundLog* checkpoint_spr_round(int invocation, int round) {
  if (invocation < 0) return nullptr;
  const CheckpointState& s = state();
  const auto it = s.spr.find({invocation, round});
  return it == s.spr.end() ? nullptr : &it->second;
}

void checkpoint_record_spr_outcome(int invocation, int round,
                                   const SprOutcome& outcome) {
  if (invocation < 0) return;
  std::ostringstream out;
  out << "P " << invocation << ' ' << round << ' ' << outcome.index << ' '
      << (outcome.moved ? 1 : 0) << ' ' << outcome.edge_a << ' '
      << outcome.edge_b << ' ' << outcome.joins[0] << ' ' << outcome.joins[1]
      << ' ' << outcome.joins[2] << ' ' << outcome.joins[3];
  record_line(out.str());
}

void checkpoint_record_spr_round_end(int invocation, int round) {
  if (invocation < 0) return;
  std::ostringstream out;
  out << "Q " << invocation << ' ' << round;
  record_line(out.str());
}

void checkpoint_note_replayed() {
  state().replayed.fetch_add(1, std::memory_order_relaxed);
}

}  // namespace qraft::internal
