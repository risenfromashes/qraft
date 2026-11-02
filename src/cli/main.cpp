#include "qraft/qraft.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <new>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#if defined(_WIN32)
#include <io.h>
#else
#include <sys/resource.h>
#include <unistd.h>
#endif

#include "cli_options.hpp"
#include "gene_input.hpp"
#include "qraft_version.hpp"
#include "selftest_data.hpp"

namespace {

using namespace qraft_cli;

bool stderr_is_tty() {
#if defined(_WIN32)
  return _isatty(_fileno(stderr)) != 0;
#else
  return isatty(STDERR_FILENO) != 0;
#endif
}

bool color_enabled() {
  static const bool enabled = [] {
    if (!stderr_is_tty() || std::getenv("NO_COLOR") != nullptr) return false;
    const char* term = std::getenv("TERM");
    return term == nullptr || std::string_view(term) != "dumb";
  }();
  return enabled;
}

std::string paint(const std::string& text, const char* code) {
  return color_enabled() ? "\x1b[" + std::string(code) + "m" + text + "\x1b[0m"
                         : text;
}

std::string version_line() {
  return std::string("qraft ") + QRAFT_VERSION + " (" + QRAFT_GIT_COMMIT + ")";
}

std::string command_list() {
  std::string out;
  for (const CommandSpec& c : command_table()) {
    std::string name = c.name;
    out += "  " + name + std::string(name.size() < 12 ? 12 - name.size() : 1, ' ') +
           c.summary + "\n";
  }
  return out;
}

constexpr const char* kExamples = R"(
Examples:
  qraft -i genes.tre -o species.tre
  qraft -i genes.tre -o species.tre -t 8 --mode extended
  qraft -i genes.tre -o species.tre --guide-tree backbone.tre
  qraft -i genes.tre -o species.tre --support localpp
  qraft support -i genes.tre --tree species.tre -o annotated.tre
  qraft check genes.tre
)";

std::string short_help() {
  std::string s = std::string("QRAFT ") + QRAFT_VERSION +
                  ": species trees from gene trees by quartet regrafting\n\n"
                  "Usage:\n"
                  "  qraft -i GENES [-o OUT] [options]    infer a species tree\n"
                  "  qraft COMMAND [options]\n\nCommands:\n" +
                  command_list() + "\nCommon options (every option: qraft help options):\n" +
                  render_options(kInfer, true) + kExamples +
                  "\nRun 'qraft help options' for every option, 'qraft help "
                  "TOPIC' for modes, guide,\nsupport, ablations, test or output, "
                  "and 'qraft COMMAND --help' for a command.\n"
                  "Quartet scores of given trees: use tqDist.\n";
  return s;
}

std::string command_help(const CommandSpec& c) {
  std::string s = std::string("Usage: ") + c.usage + "\n\n";
  for (const auto& line : [&] {
         std::vector<std::string> lines;
         std::stringstream in(c.description);
         std::string word;
         std::string line;
         while (in >> word) {
           if (!line.empty() && line.size() + 1 + word.size() > 78) {
             lines.push_back(line);
             line.clear();
           }
           line += (line.empty() ? "" : " ") + word;
         }
         if (!line.empty()) lines.push_back(line);
         return lines;
       }()) {
    s += line + "\n";
  }
  if (c.mask != 0) s += render_options(c.mask, false);
  return s;
}

std::string topic_help(const std::string& topic) {
  if (topic == "options") return command_help(*find_command("infer"));
  if (topic == "modes") {
    return "Modes:\n"
           "  default    MRE consensus, recursive divide with the bipartition "
           "test,\n             contextual regrafting and gated NNI "
           "(the QRAFT method)\n"
           "  extended   the same, then unbounded batched NNI and SPR "
           "refinement;\n             slower, for noisy biological data "
           "(--mode extended or --extended)\n"
           "Both give the same tree for any number of threads.\n" +
           render_options(kInfer, false, "Run") +
           render_options(kInfer, false, "Extended mode");
  }
  if (topic == "guide") {
    return "A guide tree (for example a trusted backbone on some taxa) edits "
           "the starting\ntree before the run.\n" +
           render_options(kInfer, false, "Guide tree");
  }
  if (topic == "support") {
    return "Local posterior probabilities can be written on the inferred tree "
           "(--support\nlocalpp) or on a given tree (qraft support).\n" +
           render_options(kInfer, false, "Branch support") + "\n" +
           command_help(*find_command("support"));
  }
  if (topic == "ablations") return render_options(kInfer, false, "Ablations");
  if (topic == "test") return render_options(kInfer, false, "Bipartition test");
  if (topic == "output") {
    return "Output:\n"
           "  The species tree is written as one Newick line to --out or to "
           "standard output.\n"
           "  --stats-json writes the options, the version and the run "
           "statistics as JSON.\n"
           "  Progress and the run summary go to standard error, coloured only "
           "on a\n  terminal and never when NO_COLOR is set; -q prints errors "
           "only.\n"
           "  Quartet scores of given trees: use tqDist.\n\n"
           "Exit status:\n"
           "  0  success\n"
           "  1  an input cannot be read or used, or the run failed\n"
           "  2  command-line error\n"
           "  3  out of memory\n";
  }
  if (topic == "man") return man_page(QRAFT_VERSION);
  if (const CommandSpec* c = find_command(topic)) return command_help(*c);
  std::vector<std::string> names = {"options", "modes", "guide", "support",
                                    "ablations", "test", "output"};
  for (const CommandSpec& c : command_table()) names.push_back(c.name);
  const std::string close = suggest(topic, names);
  throw UsageError("unknown help topic: " + topic +
                   (close.empty() ? "" : " (did you mean " + close + "?)"));
}

int hardware_threads() {
  return std::max(1, static_cast<int>(std::thread::hardware_concurrency()));
}

double peak_memory_mb() {
#if defined(_WIN32)
  return 0.0;
#else
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0) return 0.0;
  return static_cast<double>(usage.ru_maxrss) / 1024.0;  // KiB on Linux
#endif
}

std::string mode_name(const qraft::Options& o, bool support_command) {
  if (support_command) return "support on the given tree";
  const qraft::Options d;
  const bool extended = o.nni_batch && !o.nni_move_cap && !o.nni_support_gate &&
                        o.spr_refine;
  const bool custom_extended = o.nni_batch != d.nni_batch ||
                               o.nni_move_cap != d.nni_move_cap ||
                               o.nni_support_gate != d.nni_support_gate ||
                               o.spr_refine != d.spr_refine;
  const bool ablation =
      o.divide != d.divide || o.divide_score != d.divide_score ||
      o.bipartition_test.enabled != d.bipartition_test.enabled ||
      o.regraft_score != d.regraft_score || o.nni != d.nni ||
      o.bipartition_test.coefficient != d.bipartition_test.coefficient ||
      o.bipartition_test.crossover != d.bipartition_test.crossover ||
      o.bipartition_test.proposals != d.bipartition_test.proposals ||
      o.bipartition_test.anchor_pairs != d.bipartition_test.anchor_pairs ||
      o.spr_reads != d.spr_reads || o.spr_min_clade != d.spr_min_clade ||
      o.spr_candidates != d.spr_candidates || o.spr_rounds != d.spr_rounds;
  std::string name = extended && !ablation
                         ? "extended mode"
                         : (custom_extended || ablation ? "custom settings"
                                                        : "default mode");
  if (!o.starting_tree_path.empty()) name += ", given starting tree";
  if (!o.guide_tree_path.empty()) name += ", guide tree";
  if (o.support == qraft::SupportMode::LocalPP) name += ", local PP support";
  return name;
}

void write_text(const std::string& path, const std::string& text) {
  if (path.empty()) {
    std::cout << text;
    std::cout.flush();
    if (!std::cout) throw std::runtime_error("cannot write to standard output");
    return;
  }
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("cannot write the output file '" + path + "'");
  out << text;
  out.close();
  if (!out) throw std::runtime_error("cannot write the output file '" + path + "'");
}

// infer and support.
int run_inference(int begin, int argc, char** argv, bool support_command) {
  const CommandSpec& spec = *find_command(support_command ? "support" : "infer");
  CliState s;
  s.options.progress = qraft::ProgressMode::Auto;
  parse_options(begin, argc, argv, support_command ? kSupport : kInfer, s);
  if (s.help) {
    std::cout << command_help(spec);
    return 0;
  }
  if (s.genes_path.empty() && !s.positionals.empty()) {
    s.genes_path = s.positionals.front();
    s.positionals.erase(s.positionals.begin());
  }
  if (!s.positionals.empty()) {
    throw UsageError("unexpected argument: " + s.positionals.front());
  }
  qraft::Options& options = s.options;
  if (s.genes_path.empty()) throw UsageError("--genes is required");
  if (support_command) {
    if (s.tree_path.empty()) {
      throw UsageError("--tree is required: the species tree to annotate");
    }
    options.starting_tree_path = s.tree_path;
    options.divide = false;
    options.nni = false;
    options.support = qraft::SupportMode::LocalPP;
  }
  if (options.support_detail && options.support == qraft::SupportMode::None) {
    throw UsageError("--support-detail needs --support localpp");
  }
  if (options.guide_seed &&
      (options.guide_tree_path.empty() || !options.starting_tree_path.empty())) {
    throw UsageError(
        "--guide-seed needs --guide-tree and builds the MRE, so not with "
        "--starting-tree");
  }
  if (s.verbosity == Verbosity::Quiet) {
    options.progress = qraft::ProgressMode::Off;
  } else if (!s.progress_given) {
    options.progress =
        stderr_is_tty() ? qraft::ProgressMode::On : qraft::ProgressMode::Off;
  }
  // Option values the library would reject are usage errors, reported
  // before any input is read.
  try {
    qraft::validate_options(options);
  } catch (const std::invalid_argument& ex) {
    throw UsageError(ex.what());
  }
  // Every file is checked before the run, so a long run never ends on a
  // path typo.
  require_readable(s.genes_path, "gene trees file", true);
  require_regular_gene_file(s.genes_path);
  if (!options.starting_tree_path.empty()) {
    require_readable(options.starting_tree_path,
                     support_command ? "species tree file" : "starting tree file",
                     true);
  }
  if (!options.guide_tree_path.empty()) {
    require_readable(options.guide_tree_path, "guide tree file", true);
  }
  if (!options.starting_tree_path.empty()) {
    check_tree_file(options.starting_tree_path,
                    support_command ? "species tree file" : "starting tree file");
  }
  if (!options.guide_tree_path.empty()) {
    check_tree_file(options.guide_tree_path, "guide tree file");
  }
  if (!s.out_path.empty()) require_writable(s.out_path, "output file");
  if (!options.stats_json.empty()) {
    require_writable(options.stats_json, "statistics file");
  }
  if (!options.write_starting_tree_path.empty()) {
    require_writable(options.write_starting_tree_path, "starting tree output file");
  }
  options.program_version = std::string(QRAFT_VERSION) + " (" + QRAFT_GIT_COMMIT + ")";

  const int threads = options.threads > 0 ? options.threads : hardware_threads();
  const bool summary = s.verbosity == Verbosity::Verbose ||
                       (s.verbosity == Verbosity::Normal && stderr_is_tty());
  if (s.verbosity == Verbosity::Verbose) {
    std::cerr << paint("qraft:", "1") << " version " << QRAFT_VERSION << " (" << QRAFT_GIT_COMMIT << "), reading '"
              << s.genes_path << "' on " << threads << " thread"
              << (threads == 1 ? "" : "s") << ", "
              << mode_name(options, support_command) << "\n";
  }
  const auto started = std::chrono::steady_clock::now();
  std::vector<qraft::GeneTree> genes = read_gene_trees(s.genes_path, threads);
  const std::size_t gene_count = genes.size();
  const qraft::SpeciesTree tree =
      qraft::infer_species_tree(std::move(genes), options);
  write_text(s.out_path, qraft::to_newick(tree) + "\n");
  if (summary) {
    const double seconds = std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - started)
                               .count();
    std::ostringstream line;
    line << std::fixed << std::setprecision(1) << "done in " << seconds
         << " s: " << tree.leaf_labels().size() << " taxa, " << gene_count
         << " gene trees, " << mode_name(options, support_command) << ", "
         << threads << " thread" << (threads == 1 ? "" : "s")
         << ", peak memory " << std::setprecision(0) << peak_memory_mb()
         << " MB; tree written to "
         << (s.out_path.empty() ? "standard output" : "'" + s.out_path + "'");
    std::cerr << paint("qraft:", "1;32") << " " << line.str() << "\n";
  }
  return 0;
}

std::string json_string(const std::string& text) {
  std::string out = "\"";
  for (char ch : text) {
    if (ch == '"' || ch == '\\') {
      out += '\\';
      out += ch;
    } else if (static_cast<unsigned char>(ch) < 0x20) {
      char buf[8];
      std::snprintf(buf, sizeof buf, "\\u%04x", ch);
      out += buf;
    } else {
      out += ch;
    }
  }
  return out + "\"";
}

int run_check(int begin, int argc, char** argv) {
  CliState s;
  parse_options(begin, argc, argv, kCheck, s);
  if (s.help) {
    std::cout << command_help(*find_command("check"));
    return 0;
  }
  if (s.genes_path.empty() && !s.positionals.empty()) {
    s.genes_path = s.positionals.front();
    s.positionals.erase(s.positionals.begin());
  }
  if (!s.positionals.empty()) {
    throw UsageError("unexpected argument: " + s.positionals.front());
  }
  if (s.genes_path.empty()) throw UsageError("--genes is required");
  require_readable(s.genes_path, "gene trees file", true);
  require_regular_gene_file(s.genes_path);
  if (!s.json_path.empty()) require_writable(s.json_path, "summary file");
  const int threads =
      s.options.threads > 0 ? s.options.threads : hardware_threads();
  const std::vector<qraft::GeneTree> genes = read_gene_trees(s.genes_path, threads);

  std::map<std::string, std::size_t> taxon_genes;
  std::size_t min_taxa = static_cast<std::size_t>(-1);
  std::size_t max_taxa = 0;
  double sum_taxa = 0.0;
  std::size_t repeated = 0;
  std::size_t polytomies = 0;
  std::size_t small = 0;
  std::vector<std::size_t> distinct_counts;
  distinct_counts.reserve(genes.size());
  for (const qraft::GeneTree& gene : genes) {
    std::vector<std::string> labels = gene.leaf_labels();
    std::sort(labels.begin(), labels.end());
    const std::size_t leaves = labels.size();
    labels.erase(std::unique(labels.begin(), labels.end()), labels.end());
    if (labels.size() < leaves) ++repeated;
    for (const std::string& label : labels) ++taxon_genes[label];
    const std::size_t distinct = labels.size();
    distinct_counts.push_back(distinct);
    min_taxa = std::min(min_taxa, distinct);
    max_taxa = std::max(max_taxa, distinct);
    sum_taxa += static_cast<double>(distinct);
    if (distinct < 4) ++small;
    for (const qraft::Tree::Node& node : gene.nodes) {
      if (node.neighbors.size() >= 4) {
        ++polytomies;
        break;
      }
    }
  }
  const std::size_t taxa = taxon_genes.size();
  std::size_t missing = 0;
  for (std::size_t d : distinct_counts) missing += d < taxa ? 1 : 0;
  std::size_t rare = 0;
  for (const auto& [taxon, count] : taxon_genes) rare += count < 2 ? 1 : 0;
  const double mean_taxa = sum_taxa / static_cast<double>(genes.size());
  const auto percent = [&](std::size_t n) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(1)
        << 100.0 * static_cast<double>(n) / static_cast<double>(genes.size()) << "%";
    return out.str();
  };
  const bool usable = taxa >= 4 && small < genes.size();

  std::ostringstream text;
  text << s.genes_path << ": " << genes.size() << " gene trees on " << taxa
       << " taxa\n";
  text << std::fixed << std::setprecision(1) << "  taxa per gene: min "
       << min_taxa << ", mean " << mean_taxa << ", max " << max_taxa << "\n";
  text << "  genes missing taxa: " << missing << " (" << percent(missing) << ")\n";
  text << "  genes with a repeated taxon: " << repeated << " (" << percent(repeated)
       << ")\n";
  text << "  genes with polytomies: " << polytomies << " (" << percent(polytomies)
       << ")\n";
  text << "  genes with fewer than 4 taxa: " << small << " (" << percent(small)
       << ")\n";
  text << "  taxa in only one gene: " << rare << "\n";
  if (missing > 0) {
    text << "note: genes miss taxa; the bipartition test compares taxa on the "
            "genes they share.\n";
  }
  if (repeated > 0) {
    text << "note: genes repeat taxa; QRAFT uses every copy, and the "
            "bipartition test is off\n      for this input.\n";
  }
  if (small > 0) {
    text << "note: genes with fewer than 4 taxa carry no quartets and are "
            "skipped by the\n      consensus.\n";
  }
  if (usable) {
    text << "OK: 'qraft -i " << s.genes_path << "' can use this file.\n";
  }
  if (s.verbosity != Verbosity::Quiet) std::cout << text.str();
  if (!s.json_path.empty()) {
    std::ostringstream json;
    json << std::fixed << std::setprecision(3);
    json << "{\n  \"file\": " << json_string(s.genes_path) << ",\n"
         << "  \"genes\": " << genes.size() << ",\n"
         << "  \"taxa\": " << taxa << ",\n"
         << "  \"taxa_per_gene_min\": " << min_taxa << ",\n"
         << "  \"taxa_per_gene_mean\": " << mean_taxa << ",\n"
         << "  \"taxa_per_gene_max\": " << max_taxa << ",\n"
         << "  \"genes_missing_taxa\": " << missing << ",\n"
         << "  \"genes_with_repeated_taxa\": " << repeated << ",\n"
         << "  \"genes_with_polytomies\": " << polytomies << ",\n"
         << "  \"genes_below_4_taxa\": " << small << ",\n"
         << "  \"taxa_in_one_gene\": " << rare << ",\n"
         << "  \"usable\": " << (usable ? "true" : "false") << "\n}\n";
    write_text(s.json_path, json.str());
  }
  if (!usable) {
    throw std::runtime_error(
        taxa < 4 ? "the gene trees hold fewer than 4 taxa in total; QRAFT needs "
                   "at least 4"
                 : "no gene tree holds 4 or more taxa, so there are no quartets");
  }
  return 0;
}

int run_version() {
  std::cout << version_line() << "\n"
            << "commit " << QRAFT_GIT_COMMIT << ", dated " << QRAFT_BUILD_DATE
            << "\n"
            << "build: " << QRAFT_BUILD_INFO << "\n"
            << "license: AGPL-3.0-or-later\n";
  return 0;
}

std::vector<qraft::GeneTree> selftest_genes() {
  std::vector<qraft::GeneTree> genes;
  std::stringstream in(kSelftestGenes);
  for (std::string line; std::getline(in, line);) {
    if (!line.empty()) genes.push_back(qraft::parse_newick(line));
  }
  return genes;
}

int run_selftest(int begin, int argc, char** argv) {
  int threads = hardware_threads();
  for (int i = begin; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-t" || arg == "--threads") {
      if (i + 1 >= argc) throw UsageError("missing value for " + arg);
      const int given = parse_int(arg, argv[++i]);
      threads = given > 0 ? given : hardware_threads();
    } else if (arg == "-h" || arg == "--help") {
      std::cout << command_help(*find_command("selftest"));
      return 0;
    } else {
      throw UsageError("unknown argument: " + arg);
    }
  }
  struct Case {
    const char* name;
    int threads;
    bool extended;
    bool support;
    const char* expected;
  };
  const Case cases[] = {
      {"default mode", 1, false, false, kSelftestDefault},
      {"default mode", threads, false, false, kSelftestDefault},
      {"extended mode", threads, true, false, kSelftestExtended},
      {"local PP support", threads, false, true, kSelftestSupport},
  };
  int failures = 0;
  for (const Case& c : cases) {
    qraft::Options options;
    options.threads = c.threads;
    if (c.extended) {
      options.nni_batch = true;
      options.nni_move_cap = false;
      options.nni_support_gate = false;
      options.spr_refine = true;
    }
    if (c.support) options.support = qraft::SupportMode::LocalPP;
    const std::string got =
        qraft::to_newick(qraft::infer_species_tree(selftest_genes(), options)) + "\n";
    const bool ok = got == c.expected;
    failures += ok ? 0 : 1;
    std::cout << (ok ? paint("ok  ", "32") : paint("FAIL", "1;31")) << "  "
              << c.name << " (" << c.threads << " thread"
              << (c.threads == 1 ? "" : "s") << ")\n";
    if (!ok) std::cout << "      expected " << c.expected << "      got      " << got;
  }
  std::cout << (failures == 0 ? "selftest passed\n" : "selftest FAILED\n");
  return failures == 0 ? 0 : 1;
}

int run_completion(int begin, int argc, char** argv) {
  if (begin >= argc) throw UsageError("completion needs a shell: bash, zsh or fish");
  const std::string shell = argv[begin];
  if (shell == "-h" || shell == "--help") {
    std::cout << command_help(*find_command("completion"));
    return 0;
  }
  std::cout << completion_script(shell);
  return 0;
}

int dispatch(int argc, char** argv) {
  if (argc < 2) {
    std::cerr << short_help();
    return 2;
  }
  const std::string first = argv[1];
  if (first == "--version") {
    std::cout << version_line() << "\n";
    return 0;
  }
  if (first == "-h" || first == "--help") {
    std::cout << short_help();
    return 0;
  }
  if (const CommandSpec* c = find_command(first)) {
    const std::string name = c->name;
    if (name == "infer") return run_inference(2, argc, argv, false);
    if (name == "support") return run_inference(2, argc, argv, true);
    if (name == "check") return run_check(2, argc, argv);
    if (name == "version") return run_version();
    if (name == "selftest") return run_selftest(2, argc, argv);
    if (name == "completion") return run_completion(2, argc, argv);
    if (name == "help") {
      std::cout << (argc > 2 ? topic_help(argv[2]) : short_help());
      return 0;
    }
  }
  if (first[0] != '-') {
    // A first word that is no command and no file is most likely a mistyped
    // command.
    std::error_code ec;
    if (!std::filesystem::exists(first, ec)) {
      std::vector<std::string> names;
      for (const CommandSpec& c : command_table()) names.push_back(c.name);
      const std::string close = suggest(first, names);
      if (!close.empty()) {
        throw UsageError("unknown command: " + first + " (did you mean " + close +
                         "?)");
      }
    }
  }
  return run_inference(1, argc, argv, false);
}

}  // namespace

int main(int argc, char** argv) {
  try {
    return dispatch(argc, argv);
  } catch (const UsageError& ex) {
    std::cerr << paint("qraft:", "1;31") << " " << ex.what()
              << "\nRun 'qraft --help' for usage.\n";
    return 2;
  } catch (const std::bad_alloc&) {
    std::cerr << paint("qraft:", "1;31")
              << " out of memory; try fewer threads or a machine with more "
                 "memory\n";
    return 3;
  } catch (const std::exception& ex) {
    std::cerr << paint("qraft:", "1;31") << " " << ex.what() << "\n";
    return 1;
  }
}
