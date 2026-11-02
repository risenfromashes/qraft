#include "qraft/qraft.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#if defined(_WIN32)
#include <io.h>
#else
#include <unistd.h>
#endif

namespace {

constexpr const char* kUsage = R"(Usage: qraft -i GENES [-o OUT] [options]

Infers a species tree from gene trees, read as one Newick tree per line.

  -i, --genes PATH          input gene trees
  -o, --out PATH            output species tree (default: standard output)
  -t, --threads N           worker threads (default: every hardware thread)
      --starting-tree PATH  use this tree instead of the MRE consensus
      --stats-json PATH     write the options and run statistics as JSON
      --progress            report progress on standard error (default when
                            it is a terminal)
      --no-progress         do not report progress
  -h, --help                show this help

Guide tree:
      --guide-tree PATH     edit the starting tree with a guide tree on some or
                            all taxa: contract its edges that conflict with the
                            guide on the shared taxa, then refine its
                            polytomies with the guide bipartitions that keep
                            their arms whole; the run then goes on as usual
      --guide-seed          with --guide-tree and no --starting-tree, build the
                            MRE consensus without the clusters that conflict
                            with the guide, then refine it with the guide
      --write-starting-tree PATH
                            write the starting tree used, after --guide-tree

Ablations:
      --no-divide           refine the starting tree directly
      --divide-score quartet|balanced
                            rank bipartitions by satisfied-minus-violated
                            quartets (default) or pick the most balanced
      --no-bipartition-test accept the top-ranked bipartition untested
      --regraft-score raw|normalized
                            count each vertex's regrafting quartets raw
                            (default) or normalize them by its arm sizes
      --no-nni              skip gated NNI refinement
      --no-nni-cap          apply NNI moves until none qualifies (default:
                            at most 20)
      --no-nni-support-gate apply NNI moves whatever the support of the
                            bipartition they replace
      --nni-batch           apply, per scoring pass, the NNI moves of
                            largest total gain whose gains add exactly

Extended mode (unbounded batched NNI, then SPR refinement):
      --extended            --nni-batch --no-nni-cap --no-nni-support-gate
                            --spr-refine
      --spr-refine          after NNI, re-join the clades ranked by their
                            skeleton gain by the exact two-sided regraft,
                            keeping every improving move (default: off)
      --spr-reads S         core vertices read per ranking pass, more only
                            when more tie at the largest edge (default 128)
      --spr-min-clade M     fewest taxa of a ranked clade (default 10, at
                            least 2)
      --spr-candidates K    ranked clades re-joined per round (default 150)
      --spr-rounds R        rounds; a round without a move ends the stage
                            (default 2)

Branch support (on the internal edges of the output tree):
      --support localpp     write local posterior probability:length in
                            coalescent units on every internal edge whose ends
                            are binary, from one more counting pass (default:
                            none)
      --support-frequencies gene|quartet
                            quartet frequencies as ASTER counts them, each
                            gene's normalized counts summed (default;
                            ASTRAL-III also counts the unresolved quartets of
                            gene trees with polytomies), or the counts summed
                            over genes and averaged over the quartets around
                            the edge
      --support-detail      write f1-f3, EN, q1, pp1-pp3 and the length as a
                            bracketed label instead of the bare support

Bipartition test (when the first L bipartitions tried are all rejected, the
subproblem splits off one taxon):
      --test-coefficient C  threshold coefficient (default 0.05)
      --test-crossover K0   gene count at which the threshold changes sign
                            (default 10000)
      --test-proposals L    bipartitions tried per subproblem, 1-16 (default 1)
      --test-anchor-pairs H anchor pairs per bipartition, 1-4 (default 1)
)";

/// A command-line error, reported with a pointer to --help and exit status 2.
struct UsageError : std::invalid_argument {
  using std::invalid_argument::invalid_argument;
};

struct Args {
  std::string genes_path;
  std::string out_path;
  qraft::Options options;
};

bool stderr_is_tty() {
#if defined(_WIN32)
  return _isatty(_fileno(stderr)) != 0;
#else
  return isatty(STDERR_FILENO) != 0;
#endif
}

int parse_int(const std::string& flag, const std::string& value) {
  std::size_t used = 0;
  try {
    const int parsed = std::stoi(value, &used);
    if (used == value.size()) return parsed;
  } catch (const std::exception&) {
  }
  throw UsageError("invalid value for " + flag + ": " + value);
}

double parse_double(const std::string& flag, const std::string& value) {
  std::size_t used = 0;
  try {
    const double parsed = std::stod(value, &used);
    if (used == value.size()) return parsed;
  } catch (const std::exception&) {
  }
  throw UsageError("invalid value for " + flag + ": " + value);
}

Args parse_args(int argc, char** argv) {
  Args args;
  qraft::Options& options = args.options;
  options.progress = qraft::ProgressMode::Auto;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    const auto value = [&]() -> std::string {
      if (i + 1 >= argc) throw UsageError("missing value for " + arg);
      return argv[++i];
    };
    if (arg == "-i" || arg == "--genes") {
      args.genes_path = value();
    } else if (arg == "-o" || arg == "--out") {
      args.out_path = value();
    } else if (arg == "-t" || arg == "--threads") {
      options.threads = parse_int(arg, value());
    } else if (arg == "--starting-tree") {
      options.starting_tree_path = value();
    } else if (arg == "--guide-tree") {
      options.guide_tree_path = value();
    } else if (arg == "--guide-seed") {
      options.guide_seed = true;
    } else if (arg == "--write-starting-tree") {
      options.write_starting_tree_path = value();
    } else if (arg == "--stats-json") {
      options.stats_json = value();
    } else if (arg == "--progress") {
      options.progress = qraft::ProgressMode::On;
    } else if (arg == "--no-progress") {
      options.progress = qraft::ProgressMode::Off;
    } else if (arg == "--no-divide") {
      options.divide = false;
    } else if (arg == "--divide-score") {
      const std::string score = value();
      if (score == "quartet") {
        options.divide_score = qraft::DivideScore::Quartet;
      } else if (score == "balanced") {
        options.divide_score = qraft::DivideScore::Balanced;
      } else {
        throw UsageError("invalid value for --divide-score: " + score);
      }
    } else if (arg == "--no-bipartition-test") {
      options.bipartition_test.enabled = false;
    } else if (arg == "--regraft-score") {
      const std::string score = value();
      if (score == "normalized") {
        options.regraft_score = qraft::RegraftScore::Normalized;
      } else if (score == "raw") {
        options.regraft_score = qraft::RegraftScore::Raw;
      } else {
        throw UsageError("invalid value for --regraft-score: " + score);
      }
    } else if (arg == "--no-nni") {
      options.nni = false;
    } else if (arg == "--no-nni-cap") {
      options.nni_move_cap = false;
    } else if (arg == "--no-nni-support-gate") {
      options.nni_support_gate = false;
    } else if (arg == "--nni-batch") {
      options.nni_batch = true;
    } else if (arg == "--extended") {
      options.nni_batch = true;
      options.nni_move_cap = false;
      options.nni_support_gate = false;
      options.spr_refine = true;
    } else if (arg == "--spr-refine") {
      options.spr_refine = true;
    } else if (arg == "--spr-reads") {
      options.spr_reads = parse_int(arg, value());
    } else if (arg == "--spr-min-clade") {
      options.spr_min_clade = parse_int(arg, value());
    } else if (arg == "--spr-candidates") {
      options.spr_candidates = parse_int(arg, value());
    } else if (arg == "--spr-rounds") {
      options.spr_rounds = parse_int(arg, value());
    } else if (arg == "--support") {
      const std::string mode = value();
      if (mode == "localpp") {
        options.support = qraft::SupportMode::LocalPP;
      } else if (mode == "none") {
        options.support = qraft::SupportMode::None;
      } else {
        throw UsageError("invalid value for --support: " + mode);
      }
    } else if (arg == "--support-frequencies") {
      const std::string kind = value();
      if (kind == "quartet") {
        options.support_frequencies = qraft::SupportFrequencies::Quartet;
      } else if (kind == "gene") {
        options.support_frequencies = qraft::SupportFrequencies::Gene;
      } else {
        throw UsageError("invalid value for --support-frequencies: " + kind);
      }
    } else if (arg == "--support-detail") {
      options.support_detail = true;
    } else if (arg == "--test-coefficient") {
      options.bipartition_test.coefficient = parse_double(arg, value());
    } else if (arg == "--test-crossover") {
      options.bipartition_test.crossover = parse_double(arg, value());
    } else if (arg == "--test-proposals") {
      options.bipartition_test.proposals = parse_int(arg, value());
    } else if (arg == "--test-anchor-pairs") {
      options.bipartition_test.anchor_pairs = parse_int(arg, value());
    } else if (arg == "-h" || arg == "--help") {
      std::cout << kUsage;
      std::exit(0);
    } else {
      throw UsageError("unknown argument: " + arg);
    }
  }
  if (args.genes_path.empty()) throw UsageError("--genes is required");
  if (options.support_detail && options.support == qraft::SupportMode::None) {
    throw UsageError("--support-detail needs --support localpp");
  }
  if (options.guide_seed && (options.guide_tree_path.empty() || !options.starting_tree_path.empty())) {
    throw UsageError("--guide-seed needs --guide-tree and builds the MRE, so not with --starting-tree");
  }
  if (options.progress == qraft::ProgressMode::Auto) {
    options.progress =
        stderr_is_tty() ? qraft::ProgressMode::On : qraft::ProgressMode::Off;
  }
  return args;
}

bool blank(std::string_view line) {
  return std::all_of(line.begin(), line.end(), [](char ch) {
    return std::isspace(static_cast<unsigned char>(ch)) != 0;
  });
}

// The genes file is read and parsed in blocks of this many bytes, so the text
// is never held in full next to the parsed trees.
constexpr std::size_t kGeneBlockBytes = std::size_t{256} << 20;
// A block is read by one thread per this many bytes, up to `threads`.
constexpr std::size_t kGeneReadBytesPerThread = std::size_t{8} << 20;

// Runs work(t) for t in [0, count) on `count` threads, the calling thread
// being thread 0.
template <class F>
void run_on_threads(int count, F&& work) {
  if (count <= 1) {
    work(0);
    return;
  }
  // jthreads join on unwinding, so a failed launch cannot outlive the work.
  std::vector<std::jthread> pool;
  pool.reserve(static_cast<std::size_t>(count - 1));
  for (int t = 1; t < count; ++t) pool.emplace_back([&work, t] { work(t); });
  work(0);
}

// Reads bytes [offset, offset + size) of `path` into `out` on up to
// `threads` threads, each with its own stream.
void read_file_range(const std::string& path,
                     std::uint64_t offset,
                     char* out,
                     std::size_t size,
                     int threads) {
  const int readers = static_cast<int>(std::clamp<std::size_t>(
      size / kGeneReadBytesPerThread, 1, static_cast<std::size_t>(threads)));
  std::vector<char> ok(static_cast<std::size_t>(readers), 0);
  run_on_threads(readers, [&](int t) {
    const std::size_t first = size * static_cast<std::size_t>(t) /
                              static_cast<std::size_t>(readers);
    const std::size_t last = size * static_cast<std::size_t>(t + 1) /
                             static_cast<std::size_t>(readers);
    std::ifstream in(path, std::ios::binary);
    if (!in) return;
    in.seekg(static_cast<std::streamoff>(offset + first));
    in.read(out + first, static_cast<std::streamsize>(last - first));
    ok[static_cast<std::size_t>(t)] =
        in.gcount() == static_cast<std::streamsize>(last - first) ? 1 : 0;
  });
  for (char success : ok) {
    if (success == 0) throw std::runtime_error("failed to read genes file: " + path);
  }
}

// The non-blank lines of text[0, size), which ends at a line end, as views
// in order; found on up to `threads` threads.
std::vector<std::string_view> split_gene_lines(const char* text,
                                               std::size_t size,
                                               int threads) {
  const int parts = static_cast<int>(std::clamp<std::size_t>(
      size / kGeneReadBytesPerThread, 1, static_cast<std::size_t>(threads)));
  // Part t holds the lines that start in its byte range.
  std::vector<std::vector<std::string_view>> found(
      static_cast<std::size_t>(parts));
  run_on_threads(parts, [&](int t) {
    std::size_t pos = size * static_cast<std::size_t>(t) /
                      static_cast<std::size_t>(parts);
    const std::size_t end = size * static_cast<std::size_t>(t + 1) /
                            static_cast<std::size_t>(parts);
    // A line starts at 0 or after a line end.
    if (t > 0) {
      while (pos < end && text[pos - 1] != '\n') ++pos;
    }
    while (pos < end) {
      const void* nl = std::memchr(text + pos, '\n', size - pos);
      const std::size_t stop =
          nl == nullptr ? size
                        : static_cast<std::size_t>(
                              static_cast<const char*>(nl) - text);
      const std::string_view line(text + pos, stop - pos);
      if (!blank(line)) found[static_cast<std::size_t>(t)].push_back(line);
      pos = stop + 1;
    }
  });
  std::size_t count = 0;
  for (const auto& part : found) count += part.size();
  std::vector<std::string_view> lines;
  lines.reserve(count);
  for (const auto& part : found) lines.insert(lines.end(), part.begin(), part.end());
  return lines;
}

// Parses `lines` into `genes`, appended in order, on up to `threads`
// threads. The error of the first failing line is rethrown, whichever
// thread saw it.
void parse_gene_lines(const std::vector<std::string_view>& lines,
                      int threads,
                      std::vector<qraft::GeneTree>& genes) {
  const std::size_t base = genes.size();
  genes.resize(base + lines.size());
  std::vector<std::exception_ptr> failures(lines.size());
  std::atomic<std::size_t> next{0};
  const int workers = static_cast<int>(std::clamp<std::size_t>(
      lines.size(), 1, static_cast<std::size_t>(threads)));
  run_on_threads(workers, [&](int) {
    for (std::size_t i; (i = next.fetch_add(1, std::memory_order_relaxed)) <
                        lines.size();) {
      try {
        genes[base + i] = qraft::parse_newick(lines[i]);
      } catch (...) {
        failures[i] = std::current_exception();
      }
    }
  });
  for (const std::exception_ptr& failure : failures) {
    if (failure != nullptr) std::rethrow_exception(failure);
  }
}

// Parses the non-blank lines of `path`, one gene tree per line, in input
// order on `threads` threads. The file is read in blocks, each read, split
// into lines and parsed in parallel; a line cut by a block end is carried
// into the next block.
std::vector<qraft::GeneTree> read_gene_trees(const std::string& path,
                                             int threads) {
  std::ifstream probe(path, std::ios::binary | std::ios::ate);
  if (!probe) throw std::runtime_error("failed to open genes file: " + path);
  const std::streamoff file_size = probe.tellg();
  probe.close();
  if (file_size < 0) throw std::runtime_error("failed to read genes file: " + path);
  threads = std::max(1, threads);

  std::vector<qraft::GeneTree> genes;
  std::string block;
  std::size_t carried = 0;
  std::uint64_t offset = 0;
  const std::uint64_t total = static_cast<std::uint64_t>(file_size);
  while (true) {
    const std::size_t fresh = static_cast<std::size_t>(std::min<std::uint64_t>(
        kGeneBlockBytes, total - offset));
    // The carried partial line is at the front of the block.
    block.resize(carried + fresh);
    if (fresh > 0) {
      read_file_range(path, offset, block.data() + carried, fresh, threads);
    }
    offset += fresh;
    // Lines up to the last line end are complete; at the end of the file the
    // rest is the last line.
    const bool last = offset >= total;
    std::size_t complete = block.size();
    if (!last) {
      const std::size_t nl = block.rfind('\n');
      complete = nl == std::string::npos ? 0 : nl + 1;
    }
    if (complete > 0) {
      parse_gene_lines(split_gene_lines(block.data(), complete, threads),
                       threads, genes);
    }
    if (last) break;
    carried = block.size() - complete;
    std::memmove(block.data(), block.data() + complete, carried);
  }
  if (genes.empty()) {
    throw std::runtime_error("genes file did not contain any Newick trees");
  }
  if (genes.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("gene count exceeds the integer-indexed domain");
  }
  return genes;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const Args args = parse_args(argc, argv);
    const int threads =
        args.options.threads > 0
            ? args.options.threads
            : static_cast<int>(std::thread::hardware_concurrency());
    std::vector<qraft::GeneTree> genes =
        read_gene_trees(args.genes_path, threads);
    const std::string output =
        qraft::to_newick(
            qraft::infer_species_tree(std::move(genes), args.options)) +
        "\n";
    if (args.out_path.empty()) {
      std::cout << output;
    } else {
      std::ofstream out(args.out_path);
      if (!out) {
        throw std::runtime_error("failed to open output file: " + args.out_path);
      }
      out << output;
    }
    return 0;
  } catch (const UsageError& ex) {
    std::cerr << "qraft: " << ex.what() << "\nRun 'qraft --help' for usage.\n";
    return 2;
  } catch (const std::exception& ex) {
    std::cerr << "qraft: " << ex.what() << "\n";
    return 1;
  }
}
