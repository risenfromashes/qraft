#include "internal.hpp"

#include <atomic>
#include <cctype>
#include <cmath>
#include <condition_variable>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>

#if !defined(_WIN32)
#include <sys/ioctl.h>
#include <unistd.h>
#endif
#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif

namespace qraft::internal {
namespace {

constexpr long long kMetricScale = 1000000;
static_assert(std::atomic<bool>::is_always_lock_free,
              "progress notification gate must be lock-free");

/// A live subproblem; the progress line shows the one running longest.
struct ActiveSubproblemProgress {
  int depth = 0;
  int active_taxa = 0;
  double work_weight = 0.0;
  std::string phase = "start";
  std::chrono::steady_clock::time_point started;
  std::chrono::steady_clock::time_point phase_started;
};

/// Process-wide progress state. Runs that report progress are exclusive, so
/// one run owns it at a time.
struct ProgressState {
  std::mutex mutex;
  std::atomic<bool> enabled{false};
  bool rendered = false;
  bool finished = false;
  bool tree_build_enabled = true;
  bool nni_enabled = true;
  bool solve_started = false;
  int total_splits = 0;
  std::string view = "startup";  // "startup", "mre" or "solve"

  std::string mre_stage;
  int mre_total_genes = 0;
  long long mre_done_genes = 0;
  long long mre_total_entries = 0;
  long long mre_checked_entries = 0;
  int mre_accepted = 0;
  int mre_max_accepted = 0;

  double solve_total_work = 0.0;
  double solve_done_work = 0.0;
  double solve_display_fraction = 0.0;
  int live_regrafts = 0;
  int last_regraft_left_taxa = 0;
  int last_regraft_right_taxa = 0;
  int last_split_left_size = 0;
  int last_split_right_size = 0;
  long long last_split_balance_ppm = 0;
  long long worst_split_balance_ppm = 0;
  int last_active_taxa = 0;
  std::string global_stage;
  std::string nni_mode;
  int nni_pass = 0;
  int nni_pass_cap = 0;
  bool nni_done = false;
  long long nni_edges_scored = 0;
  long long nni_moves = 0;
  double nni_gain = 0.0;
  double nni_min_gain = 0.0;
  std::map<uint64_t, ActiveSubproblemProgress> active;
  bool stop_heartbeat = false;
  // Wakes the heartbeat as soon as it is stopped, so a phase change never
  // waits out the remainder of a refresh interval.
  std::condition_variable heartbeat_wake;
  std::thread heartbeat;
  std::size_t last_width = 0;
  std::chrono::steady_clock::time_point started;
  std::chrono::steady_clock::time_point last_render;
};

ProgressState& progress_state() {
  static ProgressState state;
  return state;
}

/// Tests the lock-free notification gate before touching progress state.
bool progress_enabled(const ProgressState& state) {
  return state.enabled.load(std::memory_order_acquire);
}

/// Scales a ratio to an integer in millionths.
long long scaled_metric(double value) {
  if (!std::isfinite(value) || value <= 0.0) return 0;
  return static_cast<long long>(std::llround(value * kMetricScale));
}

/// The bipartition balance, larger side / smaller side, in millionths; 0
/// when a side is empty.
long long balance_ratio_ppm(int left_size, int right_size) {
  const int smaller = std::min(left_size, right_size);
  const int larger = std::max(left_size, right_size);
  if (smaller <= 0) return 0;
  return scaled_metric(static_cast<double>(larger) /
                       static_cast<double>(smaller));
}

/// Formats a scaled diagnostic with two fractional digits.
std::string scaled_string(long long scaled) {
  std::ostringstream out;
  out << std::fixed << std::setprecision(2)
      << static_cast<double>(scaled) / static_cast<double>(kMetricScale);
  return out.str();
}

std::string compact_progress_value(const std::string& value,
                                   std::size_t max_size) {
  if (value.size() <= max_size) return value;
  if (max_size <= 1) return value.substr(0, max_size);
  return value.substr(0, max_size - 1) + "~";
}

std::string fixed_percent(double value) {
  std::ostringstream out;
  out << std::fixed << std::setprecision(1) << value;
  return out.str();
}

std::string compact_count(long long value) {
  const long long abs_value = std::llabs(value);
  std::ostringstream out;
  if (abs_value >= 1000000000ll) {
    out << std::fixed << std::setprecision(1)
        << static_cast<double>(value) / 1000000000.0 << "G";
  } else if (abs_value >= 1000000ll) {
    out << std::fixed << std::setprecision(1)
        << static_cast<double>(value) / 1000000.0 << "M";
  } else if (abs_value >= 1000ll) {
    out << std::fixed << std::setprecision(1)
        << static_cast<double>(value) / 1000.0 << "k";
  } else {
    out << value;
  }
  return out.str();
}

std::string make_bar(double fraction, int width) {
  fraction = std::clamp(fraction, 0.0, 1.0);
  const int filled =
      std::clamp(static_cast<int>(std::llround(fraction * width)), 0, width);
  std::string bar;
  bar.reserve(static_cast<std::size_t>(width) + 2);
  bar.push_back('[');
  for (int i = 0; i < width; ++i) {
    bar.push_back(i < filled ? '#' : '-');
  }
  bar.push_back(']');
  return bar;
}

/// The progress work of a subproblem on m taxa, m^2. It is credited 32% when
/// the subproblem's bipartition is chosen and 68% when its sides are
/// regrafted.
double split_work_weight(int active_size) {
  const double m = static_cast<double>(std::max(0, active_size));
  return m * m;
}

std::string nni_pass_string(const ProgressState& state) {
  if (!state.nni_enabled) return "nni=skip";
  if (state.nni_mode.empty()) return "nni=pending";
  std::ostringstream out;
  out << "nni_pass=" << state.nni_pass;
  if (state.nni_done) out << " done";
  return out.str();
}

/// The terminal width of standard error, or 160 when unknown.
int terminal_width() {
#if !defined(_WIN32)
  winsize size{};
  if (ioctl(STDERR_FILENO, TIOCGWINSZ, &size) == 0 && size.ws_col > 0) {
    return static_cast<int>(size.ws_col);
  }
#endif
  return 160;
}

/// Reads a short CPU model string for startup diagnostics.
std::string cpu_model_name() {
#if defined(__APPLE__)
  char brand[256] = {};
  std::size_t size = sizeof brand - 1;
  if (sysctlbyname("machdep.cpu.brand_string", brand, &size, nullptr, 0) == 0 &&
      brand[0] != '\0') {
    return brand;
  }
#elif !defined(_WIN32)
  std::ifstream in("/proc/cpuinfo");
  std::string line;
  const std::string prefix = "model name";
  while (std::getline(in, line)) {
    if (line.rfind(prefix, 0) == 0) {
      const std::size_t colon = line.find(':');
      if (colon != std::string::npos) {
        std::string value = line.substr(colon + 1);
        while (!value.empty() &&
               std::isspace(static_cast<unsigned char>(value.front()))) {
          value.erase(value.begin());
        }
        return value;
      }
    }
  }
#endif
  return "unknown CPU";
}

/// Emits one carriage-return progress line and clears any previous tail.
void emit_progress_line_locked(ProgressState& state, std::string line) {
  const int width = std::max(40, terminal_width());
  if (static_cast<int>(line.size()) >= width) {
    line.resize(static_cast<std::size_t>(width - 1));
  }
  std::cerr << '\r' << line << "\x1b[K";
  std::cerr.flush();
  state.last_width = line.size();
  state.rendered = true;
}

/// Emits the progress line, at most every 100 ms unless `force`. The caller
/// holds the state mutex.
void render_progress_locked(ProgressState& state, bool force) {
  if (!progress_enabled(state)) return;
  const auto now = std::chrono::steady_clock::now();
  if (!force && now - state.last_render < std::chrono::milliseconds(100)) {
    return;
  }
  state.last_render = now;

  const double elapsed =
      std::chrono::duration<double>(now - state.started).count();
  const RunStats& stats = run_stats();
  std::ostringstream out;
  const int width = std::max(80, terminal_width());
  const int bar_width = std::clamp(width / 4, 18, 42);

  if (state.view == "mre") {
    // Split counting fills the first 65% of the bar and assembly the rest.
    double fraction = 0.0;
    if (state.mre_stage == "count") {
      const int total = std::max(1, state.mre_total_genes);
      fraction = 0.65 *
                 static_cast<double>(std::min<long long>(
                     state.mre_done_genes, total)) /
                 static_cast<double>(total);
    } else if (state.mre_stage == "assemble") {
      const double accepted_fraction =
          state.mre_max_accepted > 0
              ? static_cast<double>(std::min(state.mre_accepted,
                                             state.mre_max_accepted)) /
                    static_cast<double>(state.mre_max_accepted)
              : 1.0;
      const double checked_fraction =
          state.mre_total_entries > 0
              ? static_cast<double>(std::min<long long>(
                    state.mre_checked_entries,
                    state.mre_total_entries)) /
                    static_cast<double>(state.mre_total_entries)
              : accepted_fraction;
      fraction = 0.65 + 0.35 * std::max(accepted_fraction, checked_fraction);
    } else if (state.mre_stage == "done") {
      fraction = 1.0;
    }
    fraction = std::clamp(fraction, 0.0, 1.0);
    out << "[qraft] mre " << make_bar(fraction, bar_width) << " "
        << fixed_percent(100.0 * fraction) << "% elapsed=" << std::fixed
        << std::setprecision(1) << elapsed << "s";
    if (state.mre_stage == "count") {
      out << " count=" << state.mre_done_genes << "/"
          << state.mre_total_genes << " genes";
      // The count stage has a known total, so its remaining time can be
      // estimated from the genes done so far.
      if (state.mre_done_genes > 0 &&
          state.mre_done_genes < state.mre_total_genes) {
        const double remaining =
            elapsed *
            static_cast<double>(state.mre_total_genes - state.mre_done_genes) /
            static_cast<double>(state.mre_done_genes);
        out << " eta=" << std::setprecision(1) << remaining << "s";
      }
    } else if (state.mre_stage == "assemble") {
      out << " assemble checked="
          << compact_count(state.mre_checked_entries) << "/"
          << compact_count(state.mre_total_entries)
          << " accepted=" << state.mre_accepted << "/"
          << state.mre_max_accepted;
    } else {
      out << " done accepted=" << state.mre_accepted << "/"
          << state.mre_max_accepted;
    }
    emit_progress_line_locked(state, out.str());
    return;
  }

  const ActiveSubproblemProgress* oldest = nullptr;
  double oldest_elapsed = 0.0;
  for (const auto& entry : state.active) {
    const double active_elapsed =
        std::chrono::duration<double>(now - entry.second.started).count();
    if (oldest == nullptr || active_elapsed > oldest_elapsed) {
      oldest = &entry.second;
      oldest_elapsed = active_elapsed;
    }
  }

  double fraction = state.solve_total_work > 0.0
                        ? state.solve_done_work / state.solve_total_work
                        : (state.finished ? 1.0 : 0.0);
  if (state.finished) fraction = 1.0;
  fraction = std::clamp(fraction, 0.0, 1.0);
  // New subproblems add to the total work; the bar never moves back.
  state.solve_display_fraction =
      std::max(state.solve_display_fraction, fraction);
  const double display_fraction =
      state.finished ? 1.0 : state.solve_display_fraction;
  out << "[qraft] solve " << make_bar(display_fraction, bar_width) << " "
      << fixed_percent(100.0 * display_fraction) << "% elapsed="
      << std::fixed << std::setprecision(1) << elapsed << "s "
      << "work=" << compact_count(static_cast<long long>(
             std::llround(state.solve_done_work)))
      << "/" << compact_count(static_cast<long long>(
             std::llround(std::max(state.solve_total_work, 1.0))))
      << " ";
  if (oldest != nullptr) {
    const double phase_elapsed =
        std::chrono::duration<double>(now - oldest->phase_started).count();
    out << "live=" << state.active.size() << " d=" << oldest->depth
        << " taxa=" << oldest->active_taxa
        << " phase=" << compact_progress_value(oldest->phase, 18)
        << " " << std::setprecision(1) << phase_elapsed << "s ";
  } else {
    out << "live=0 last_taxa=" << state.last_active_taxa << " ";
  }
  if (state.last_split_left_size > 0 || state.last_split_right_size > 0) {
    out << "L=" << state.last_split_left_size
        << " R=" << state.last_split_right_size
        << " bal=" << scaled_string(state.last_split_balance_ppm) << " ";
  }
  if (state.live_regrafts > 0) {
    out << "regrafts=" << state.live_regrafts
        << " L=" << state.last_regraft_left_taxa
        << " R=" << state.last_regraft_right_taxa << " ";
  }
  if (!state.global_stage.empty()) {
    out << "stage=" << compact_progress_value(state.global_stage, 36) << " ";
  }
  out << nni_pass_string(state) << " hdt_q="
      << compact_count(stats.hdt_queries.load(std::memory_order_relaxed))
      << " worst_bal="
      << scaled_string(state.worst_split_balance_ppm);

  emit_progress_line_locked(state, out.str());
}

/// Redraws the line every second, so elapsed times stay current while a long
/// subproblem runs.
void progress_heartbeat_loop() {
  ProgressState& state = progress_state();
  std::unique_lock<std::mutex> lock(state.mutex);
  while (true) {
    state.heartbeat_wake.wait_for(lock, std::chrono::seconds(1),
                                  [&] { return state.stop_heartbeat; });
    if (state.stop_heartbeat || !progress_enabled(state)) return;
    render_progress_locked(state, true);
  }
}

/// Stops and joins the progress heartbeat thread if one is active.
void stop_progress_heartbeat() {
  ProgressState& state = progress_state();
  std::thread heartbeat;
  {
    std::lock_guard<std::mutex> lock(state.mutex);
    state.stop_heartbeat = true;
    heartbeat = std::move(state.heartbeat);
  }
  state.heartbeat_wake.notify_all();
  if (heartbeat.joinable()) heartbeat.join();
}

}  // namespace

void progress_startup_status(const Options& options,
                             std::string stage,
                             std::string detail) {
  if (options.progress != ProgressMode::On) return;
  stop_progress_heartbeat();
  ProgressState& state = progress_state();
  std::lock_guard<std::mutex> lock(state.mutex);
  if (!state.rendered) {
    state.started = std::chrono::steady_clock::now();
    state.last_render = state.started;
    state.last_width = 0;
  }
  state.enabled.store(true, std::memory_order_release);
  state.view = "startup";
  const double elapsed =
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    state.started)
          .count();
  std::ostringstream out;
  out << "[qraft] startup stage=" << stage << " elapsed=" << std::fixed
      << std::setprecision(1) << elapsed << "s"
      << " workers=" << std::max(1, get_num_workers())
      << " cpu=\"" << cpu_model_name() << "\"";
  if (!detail.empty()) {
    out << " " << detail;
  }
  emit_progress_line_locked(state, out.str());
}

/// Starts the progress bar of the MRE starting-tree build.
void progress_mre_start(const Options& options,
                              int total_genes,
                              int taxon_count) {
  if (options.progress != ProgressMode::On) return;
  (void)taxon_count;
  stop_progress_heartbeat();
  ProgressState& state = progress_state();
  bool should_start_heartbeat = false;
  {
    std::lock_guard<std::mutex> lock(state.mutex);
    state.enabled.store(options.progress == ProgressMode::On,
                        std::memory_order_release);
    state.rendered = false;
    state.finished = false;
    state.view = "mre";
    state.mre_stage = "count";
    state.mre_total_genes = std::max(0, total_genes);
    state.mre_done_genes = 0;
    state.mre_total_entries = 0;
    state.mre_checked_entries = 0;
    state.mre_accepted = 0;
    state.mre_max_accepted = std::max(0, taxon_count - 3);
    state.stop_heartbeat = false;
    state.started = std::chrono::steady_clock::now();
    state.last_render = state.started;
    render_progress_locked(state, true);
    should_start_heartbeat = progress_enabled(state);
  }
  if (should_start_heartbeat) {
    state.heartbeat = std::thread(progress_heartbeat_loop);
  }
}

/// Records one gene whose splits the MRE build has counted.
void progress_mre_count_gene() {
  ProgressState& state = progress_state();
  if (!progress_enabled(state)) return;
  std::lock_guard<std::mutex> lock(state.mutex);
  if (!progress_enabled(state) || state.view != "mre") return;
  state.mre_stage = "count";
  state.mre_done_genes =
      std::min<long long>(std::max(0, state.mre_total_genes),
                          state.mre_done_genes + 1);
  render_progress_locked(state, false);
}

/// Updates MRE tree assembly: split entries checked and splits accepted.
void progress_mre_assembly(long long checked_entries,
                                 long long total_entries,
                                 int accepted,
                                 int max_accepted) {
  ProgressState& state = progress_state();
  if (!progress_enabled(state)) return;
  std::lock_guard<std::mutex> lock(state.mutex);
  if (!progress_enabled(state) || state.view != "mre") return;
  state.mre_stage = "assemble";
  state.mre_total_entries = std::max(0ll, total_entries);
  state.mre_checked_entries = std::max(0ll, checked_entries);
  state.mre_accepted = std::max(0, accepted);
  state.mre_max_accepted = std::max(0, max_accepted);
  render_progress_locked(state, false);
}

/// Finishes the MRE progress bar.
void progress_mre_finish(int accepted, int max_accepted) {
  ProgressState& state = progress_state();
  if (!progress_enabled(state)) return;
  stop_progress_heartbeat();
  std::lock_guard<std::mutex> lock(state.mutex);
  if (!progress_enabled(state) || state.view != "mre") return;
  state.mre_stage = "done";
  state.mre_accepted = std::max(0, accepted);
  state.mre_max_accepted = std::max(0, max_accepted);
  render_progress_locked(state, true);
  if (state.rendered) {
    std::cerr << '\n';
  }
  state.enabled.store(false, std::memory_order_release);
}

/// Initializes progress reporting for one inference run.
void progress_start(const Options& options, int total_splits) {
  if (options.progress != ProgressMode::On) return;
  stop_progress_heartbeat();
  ProgressState& state = progress_state();
  bool should_start_heartbeat = false;
  {
    std::lock_guard<std::mutex> lock(state.mutex);
    state.enabled.store(options.progress == ProgressMode::On,
                        std::memory_order_release);
    state.rendered = false;
    state.finished = false;
    state.view = "solve";
    state.solve_started = true;
    state.tree_build_enabled = options.divide;
    state.nni_enabled = options.nni;
    state.total_splits = std::max(0, total_splits);
    state.mre_stage.clear();
    state.solve_total_work = 0.0;
    state.solve_done_work = 0.0;
    state.solve_display_fraction = 0.0;
    state.live_regrafts = 0;
    state.last_regraft_left_taxa = 0;
    state.last_regraft_right_taxa = 0;
    state.last_split_left_size = 0;
    state.last_split_right_size = 0;
    state.last_split_balance_ppm = 0;
    state.worst_split_balance_ppm = 0;
    state.last_active_taxa = 0;
    state.global_stage.clear();
    state.nni_mode.clear();
    state.nni_pass = 0;
    state.nni_pass_cap = 0;
    state.nni_done = false;
    state.nni_edges_scored = 0;
    state.nni_moves = 0;
    state.nni_gain = 0.0;
    state.nni_min_gain = 0.0;
    state.active.clear();
    state.stop_heartbeat = false;
    state.last_width = 0;
    state.started = std::chrono::steady_clock::now();
    state.last_render = state.started;
    render_progress_locked(state, true);
    should_start_heartbeat = progress_enabled(state);
  }
  if (should_start_heartbeat) {
    state.heartbeat = std::thread(progress_heartbeat_loop);
  }
}

/// Records one non-terminal recursive subproblem.
void progress_note_subproblem(uint64_t path, int depth, int active_taxa) {
  ProgressState& state = progress_state();
  if (!progress_enabled(state)) return;
  std::lock_guard<std::mutex> lock(state.mutex);
  if (!progress_enabled(state)) return;
  state.last_active_taxa = active_taxa;
  const auto now = std::chrono::steady_clock::now();
  const double work_weight = split_work_weight(active_taxa);
  if (state.tree_build_enabled) {
    state.solve_total_work += work_weight;
  }
  state.active[path] = ActiveSubproblemProgress{
      depth, active_taxa, work_weight, "start", now, now};
  render_progress_locked(state, false);
}

/// Updates the phase for a live recursive subproblem.
void progress_note_phase(uint64_t path, std::string phase) {
  ProgressState& state = progress_state();
  if (!progress_enabled(state)) return;
  std::lock_guard<std::mutex> lock(state.mutex);
  if (!progress_enabled(state)) return;
  const auto it = state.active.find(path);
  if (it == state.active.end()) return;
  it->second.phase = std::move(phase);
  it->second.phase_started = std::chrono::steady_clock::now();
  render_progress_locked(state, false);
}

/// Updates the global phase label.
void progress_note_stage(std::string stage) {
  ProgressState& state = progress_state();
  if (!progress_enabled(state)) return;
  std::lock_guard<std::mutex> lock(state.mutex);
  if (!progress_enabled(state)) return;
  state.global_stage = std::move(stage);
  render_progress_locked(state, false);
}

void progress_note_regraft_start(int parent_active_size,
                                 int left_taxa,
                                 int right_taxa) {
  ProgressState& state = progress_state();
  if (!progress_enabled(state)) return;
  std::lock_guard<std::mutex> lock(state.mutex);
  if (!progress_enabled(state)) return;
  state.live_regrafts = std::max(0, state.live_regrafts) + 1;
  state.last_active_taxa = std::max(0, parent_active_size);
  state.last_regraft_left_taxa = std::max(0, left_taxa);
  state.last_regraft_right_taxa = std::max(0, right_taxa);
  state.global_stage = "regraft";
  render_progress_locked(state, false);
}

void progress_note_regraft_done(int parent_active_size) {
  ProgressState& state = progress_state();
  if (!progress_enabled(state)) return;
  std::lock_guard<std::mutex> lock(state.mutex);
  if (!progress_enabled(state)) return;
  state.live_regrafts = std::max(0, state.live_regrafts - 1);
  if (state.tree_build_enabled) {
    state.solve_done_work += 0.68 * split_work_weight(parent_active_size);
  }
  render_progress_locked(state, false);
}

void progress_note_nni(std::string mode,
                           int moves,
                           int max_moves,
                           long long edges_scored,
                           double gain,
                           double min_gain) {
  ProgressState& state = progress_state();
  if (!progress_enabled(state)) return;
  std::lock_guard<std::mutex> lock(state.mutex);
  if (!progress_enabled(state)) return;
  state.nni_mode = std::move(mode);
  state.nni_pass = moves;
  state.nni_pass_cap = std::max(0, max_moves);
  state.nni_done = state.nni_mode.find("done") != std::string::npos;
  state.nni_edges_scored = std::max(0ll, edges_scored);
  state.nni_moves = std::max(0, moves);
  state.nni_gain = gain;
  state.nni_min_gain = std::max(0.0, min_gain);
  render_progress_locked(state, false);
}

void progress_note_split(uint64_t path,
                         int active_taxa,
                         int left_size,
                         int right_size) {
  ProgressState& state = progress_state();
  if (!progress_enabled(state)) return;
  std::lock_guard<std::mutex> lock(state.mutex);
  if (!progress_enabled(state)) return;
  state.last_active_taxa = active_taxa;
  const auto it = state.active.find(path);
  const double work_weight =
      it == state.active.end() ? split_work_weight(active_taxa)
                               : it->second.work_weight;
  state.active.erase(path);
  state.last_split_left_size = left_size;
  state.last_split_right_size = right_size;
  state.last_split_balance_ppm = balance_ratio_ppm(left_size, right_size);
  state.worst_split_balance_ppm =
      std::max(state.worst_split_balance_ppm,
               state.last_split_balance_ppm);
  if (state.tree_build_enabled) {
    state.solve_done_work += 0.32 * work_weight;
  }
  render_progress_locked(state, false);
}

/// Finishes progress output for one inference run.
void progress_finish() {
  ProgressState& state = progress_state();
  if (!progress_enabled(state)) return;
  stop_progress_heartbeat();
  std::lock_guard<std::mutex> lock(state.mutex);
  state.finished = true;
  state.live_regrafts = 0;
  if (state.nni_enabled) {
    state.nni_done = true;
  }
  render_progress_locked(state, true);
  if (progress_enabled(state) && state.rendered) {
    std::cerr << '\n';
  }
  state.enabled.store(false, std::memory_order_release);
  state.finished = false;
  state.active.clear();
}

}  // namespace qraft::internal
