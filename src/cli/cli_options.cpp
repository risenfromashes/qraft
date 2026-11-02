#include "cli_options.hpp"

#include <algorithm>
#include <cstddef>
#include <exception>
#include <sstream>

namespace qraft_cli {

int parse_int(const std::string& flag, const std::string& value) {
  std::size_t used = 0;
  try {
    const int parsed = std::stoi(value, &used);
    if (used == value.size()) return parsed;
  } catch (const std::exception&) {
  }
  throw UsageError("invalid value for " + flag + ": " + value +
                   " (expected a whole number)");
}

namespace {

double parse_double(const std::string& flag, const std::string& value) {
  std::size_t used = 0;
  try {
    const double parsed = std::stod(value, &used);
    if (used == value.size()) return parsed;
  } catch (const std::exception&) {
  }
  throw UsageError("invalid value for " + flag + ": " + value +
                   " (expected a number)");
}

constexpr unsigned kRun = kInfer | kSupport | kScore;
constexpr unsigned kAll = kInfer | kSupport | kScore | kCheck;

// Help sections, in order.
const char* const kGroups[] = {
    "Input and output", "Run",
    "Extended mode",    "Guide tree",
    "Branch support",   "Ablations",
    "Bipartition test",
};

const char* group_heading(const std::string& group) {
  if (group == "Extended mode") {
    return "Extended mode (unbounded batched NNI, then SPR refinement):";
  }
  if (group == "Branch support") {
    return "Branch support (on the internal edges of the output tree):";
  }
  if (group == "Bipartition test") {
    return "Bipartition test (when the first L bipartitions tried are all "
           "rejected, the\nsubproblem splits off one taxon):";
  }
  if (group == "Ablations") {
    return "Ablations (for experiments; the defaults are the QRAFT method):";
  }
  return nullptr;
}

void set_extended(qraft::Options& o, bool on) {
  o.nni_batch = on;
  o.nni_move_cap = !on;
  o.nni_support_gate = !on;
  o.spr_refine = on;
}

}  // namespace

const std::vector<OptionSpec>& option_table() {
  static const std::vector<OptionSpec> table = {
      // Input and output.
      {"--genes", "-i", "PATH", nullptr, "Input and output", kAll, true,
       "gene trees, one Newick tree per line (may also be given as the "
       "first argument)",
       [](CliState& s, const std::string& v) { s.genes_path = v; }},
      {"--out", "-o", "PATH", nullptr, "Input and output", kInfer | kSupport,
       true, "output species tree (default: standard output)",
       [](CliState& s, const std::string& v) { s.out_path = v; }},
      {"--tree", nullptr, "PATH", nullptr, "Input and output",
       kSupport | kScore, true, "the species tree to annotate or score",
       [](CliState& s, const std::string& v) { s.tree_path = v; }},
      {"--stats-json", nullptr, "PATH", nullptr, "Input and output",
       kInfer | kSupport, true,
       "write the options and run statistics as JSON",
       [](CliState& s, const std::string& v) { s.options.stats_json = v; }},
      {"--json", nullptr, "PATH", nullptr, "Input and output", kCheck, true,
       "also write the summary as JSON",
       [](CliState& s, const std::string& v) { s.json_path = v; }},
      {"--starting-tree", nullptr, "PATH", nullptr, "Input and output", kInfer,
       false, "use this tree instead of the MRE consensus",
       [](CliState& s, const std::string& v) {
         s.options.starting_tree_path = v;
       }},
      // Run.
      {"--threads", "-t", "N", nullptr, "Run", kAll, true,
       "worker threads (default 0, or any value below 1: every hardware "
       "thread); the tree does not depend on it",
       [](CliState& s, const std::string& v) {
         s.options.threads = parse_int("--threads", v);
       }},
      {"--mode", nullptr, "default|extended", "default|extended", "Run",
       kInfer, true,
       "the pipeline: default, or extended (as --extended: slower, for "
       "noisy biological data)",
       [](CliState& s, const std::string& v) {
         set_extended(s.options, v == "extended");
       }},
      {"--progress", nullptr, nullptr, nullptr, "Run", kRun, false,
       "report progress on standard error (default when it is a terminal)",
       [](CliState& s, const std::string&) {
         s.options.progress = qraft::ProgressMode::On;
         s.progress_given = true;
       }},
      {"--no-progress", nullptr, nullptr, nullptr, "Run", kRun, false,
       "do not report progress",
       [](CliState& s, const std::string&) {
         s.options.progress = qraft::ProgressMode::Off;
         s.progress_given = true;
       }},
      {"--quiet", "-q", nullptr, nullptr, "Run", kAll, true,
       "print only errors: no progress and no summary",
       [](CliState& s, const std::string&) { s.verbosity = Verbosity::Quiet; }},
      {"--verbose", "-v", nullptr, nullptr, "Run", kAll, true,
       "also print what is read and the run summary when not on a terminal",
       [](CliState& s, const std::string&) {
         s.verbosity = Verbosity::Verbose;
       }},
      {"--help", "-h", nullptr, nullptr, "Run", kAll, true,
       "show this help",
       [](CliState& s, const std::string&) { s.help = true; }},
      // Extended mode.
      {"--extended", nullptr, nullptr, nullptr, "Extended mode", kInfer, false,
       "--nni-batch --no-nni-cap --no-nni-support-gate --spr-refine",
       [](CliState& s, const std::string&) { set_extended(s.options, true); }},
      {"--spr-refine", nullptr, nullptr, nullptr, "Extended mode", kInfer,
       false,
       "after NNI, re-join the clades ranked by their skeleton gain by the "
       "exact two-sided regraft, keeping every improving move (default: off)",
       [](CliState& s, const std::string&) { s.options.spr_refine = true; }},
      {"--spr-reads", nullptr, "S", nullptr, "Extended mode", kInfer, false,
       "core vertices read per ranking pass, more only when more tie at the "
       "largest edge (default 128)",
       [](CliState& s, const std::string& v) {
         s.options.spr_reads = parse_int("--spr-reads", v);
       }},
      {"--spr-min-clade", nullptr, "M", nullptr, "Extended mode", kInfer,
       false, "fewest taxa of a ranked clade (default 10, at least 2)",
       [](CliState& s, const std::string& v) {
         s.options.spr_min_clade = parse_int("--spr-min-clade", v);
       }},
      {"--spr-candidates", nullptr, "K", nullptr, "Extended mode", kInfer,
       false, "ranked clades re-joined per round (default 150)",
       [](CliState& s, const std::string& v) {
         s.options.spr_candidates = parse_int("--spr-candidates", v);
       }},
      {"--spr-rounds", nullptr, "R", nullptr, "Extended mode", kInfer, false,
       "rounds; a round without a move ends the stage (default 2)",
       [](CliState& s, const std::string& v) {
         s.options.spr_rounds = parse_int("--spr-rounds", v);
       }},
      // Guide tree.
      {"--guide-tree", nullptr, "PATH", nullptr, "Guide tree", kInfer, true,
       "edit the starting tree with a guide tree on some or all taxa: "
       "contract its edges that conflict with the guide on the shared taxa, "
       "then refine its polytomies with the guide bipartitions that keep "
       "their arms whole; the run then goes on as usual",
       [](CliState& s, const std::string& v) { s.options.guide_tree_path = v; }},
      {"--guide-seed", nullptr, nullptr, nullptr, "Guide tree", kInfer, false,
       "with --guide-tree and no --starting-tree, build the MRE consensus "
       "without the clusters that conflict with the guide, then refine it "
       "with the guide",
       [](CliState& s, const std::string&) { s.options.guide_seed = true; }},
      {"--write-starting-tree", nullptr, "PATH", nullptr, "Guide tree", kInfer,
       false, "write the starting tree used, after --guide-tree",
       [](CliState& s, const std::string& v) {
         s.options.write_starting_tree_path = v;
       }},
      // Branch support.
      {"--support", nullptr, "localpp", "localpp|none", "Branch support",
       kInfer, true,
       "write local posterior probability:length in coalescent units on "
       "every internal edge whose ends are binary, from one more counting "
       "pass (default: none)",
       [](CliState& s, const std::string& v) {
         s.options.support = v == "localpp" ? qraft::SupportMode::LocalPP
                                            : qraft::SupportMode::None;
       }},
      {"--support-frequencies", nullptr, "gene|quartet", "gene|quartet",
       "Branch support", kInfer | kSupport, false,
       "quartet frequencies as ASTER counts them, each gene's normalized "
       "counts summed (default; ASTRAL-III also counts the unresolved "
       "quartets of gene trees with polytomies), or the counts summed over "
       "genes and averaged over the quartets around the edge",
       [](CliState& s, const std::string& v) {
         s.options.support_frequencies = v == "quartet"
                                             ? qraft::SupportFrequencies::Quartet
                                             : qraft::SupportFrequencies::Gene;
       }},
      {"--support-detail", nullptr, nullptr, nullptr, "Branch support",
       kInfer | kSupport, false,
       "write f1-f3, EN, q1, pp1-pp3 and the length as a bracketed label "
       "instead of the bare support",
       [](CliState& s, const std::string&) { s.options.support_detail = true; }},
      // Ablations.
      {"--no-divide", nullptr, nullptr, nullptr, "Ablations", kInfer, false,
       "refine the starting tree directly",
       [](CliState& s, const std::string&) { s.options.divide = false; }},
      {"--divide-score", nullptr, "quartet|balanced", "quartet|balanced",
       "Ablations", kInfer, false,
       "rank bipartitions by satisfied-minus-violated quartets (default) or "
       "pick the most balanced",
       [](CliState& s, const std::string& v) {
         s.options.divide_score = v == "balanced" ? qraft::DivideScore::Balanced
                                                  : qraft::DivideScore::Quartet;
       }},
      {"--no-bipartition-test", nullptr, nullptr, nullptr, "Ablations", kInfer,
       false, "accept the top-ranked bipartition untested",
       [](CliState& s, const std::string&) {
         s.options.bipartition_test.enabled = false;
       }},
      {"--regraft-score", nullptr, "raw|normalized", "raw|normalized",
       "Ablations", kInfer, false,
       "count each vertex's regrafting quartets raw (default) or normalize "
       "them by its arm sizes",
       [](CliState& s, const std::string& v) {
         s.options.regraft_score = v == "normalized"
                                       ? qraft::RegraftScore::Normalized
                                       : qraft::RegraftScore::Raw;
       }},
      {"--no-nni", nullptr, nullptr, nullptr, "Ablations", kInfer, false,
       "skip gated NNI refinement",
       [](CliState& s, const std::string&) { s.options.nni = false; }},
      {"--no-nni-cap", nullptr, nullptr, nullptr, "Ablations", kInfer, false,
       "apply NNI moves until none qualifies (default: at most 20)",
       [](CliState& s, const std::string&) { s.options.nni_move_cap = false; }},
      {"--no-nni-support-gate", nullptr, nullptr, nullptr, "Ablations", kInfer,
       false,
       "apply NNI moves whatever the support of the bipartition they replace",
       [](CliState& s, const std::string&) {
         s.options.nni_support_gate = false;
       }},
      {"--nni-batch", nullptr, nullptr, nullptr, "Ablations", kInfer, false,
       "apply, per scoring pass, the NNI moves of largest total gain whose "
       "gains add exactly",
       [](CliState& s, const std::string&) { s.options.nni_batch = true; }},
      {"--refine-start", nullptr, nullptr, nullptr, "Ablations", kInfer,
       false,
       "first refine the starting tree as --no-divide --extended does, then "
       "run from it (the refined MRE; the same trees as the commands in "
       "'qraft help ablations')",
       [](CliState& s, const std::string&) { s.options.refine_start = true; }},
      // Bipartition test.
      {"--test-coefficient", nullptr, "C", nullptr, "Bipartition test", kInfer,
       false, "threshold coefficient (default 0.05)",
       [](CliState& s, const std::string& v) {
         s.options.bipartition_test.coefficient =
             parse_double("--test-coefficient", v);
       }},
      {"--test-crossover", nullptr, "K0", nullptr, "Bipartition test", kInfer,
       false, "gene count at which the threshold changes sign (default 10000)",
       [](CliState& s, const std::string& v) {
         s.options.bipartition_test.crossover =
             parse_double("--test-crossover", v);
       }},
      {"--test-proposals", nullptr, "L", nullptr, "Bipartition test", kInfer,
       false, "bipartitions tried per subproblem, 1-16 (default 1)",
       [](CliState& s, const std::string& v) {
         s.options.bipartition_test.proposals =
             parse_int("--test-proposals", v);
       }},
      {"--test-anchor-pairs", nullptr, "H", nullptr, "Bipartition test",
       kInfer, false, "anchor pairs per bipartition, 1-4 (default 1)",
       [](CliState& s, const std::string& v) {
         s.options.bipartition_test.anchor_pairs =
             parse_int("--test-anchor-pairs", v);
       }},
  };
  return table;
}

const std::vector<CommandSpec>& command_table() {
  static const std::vector<CommandSpec> table = {
      {"infer", kInfer, "qraft [infer] -i GENES [-o OUT] [options]",
       "infer a species tree from gene trees (the default command)",
       "Infers a species tree from gene trees, read as one Newick tree per "
       "line. 'qraft -i GENES' and 'qraft infer -i GENES' are the same."},
      {"support", kSupport, "qraft support -i GENES --tree TREE [-o OUT]",
       "write local posterior probabilities on the edges of a given tree",
       "Reads a species tree and writes it back with the local posterior "
       "probability and the length in coalescent units of every internal edge "
       "whose ends are binary, from the quartet frequencies of the gene trees. "
       "The topology is not changed. Same as 'qraft -i GENES --starting-tree "
       "TREE --no-divide --no-nni --support localpp'."},
      {"check", kCheck, "qraft check [-i] GENES [--json PATH]",
       "check a gene-tree file and summarize it, without inferring",
       "Reads the gene trees and reports how many there are, their taxa, "
       "missing and repeated taxa, polytomies and genes too small to use, "
       "and what QRAFT will do with them. Exits with status 1 when the file "
       "cannot be used."},
      {"version", 0, "qraft version", "print the version and build details",
       "Prints the version, the git commit, the compiler and the build "
       "options. 'qraft --version' prints the first line only."},
      {"selftest", 0, "qraft selftest [-t N]",
       "run a small built-in example and check the result",
       "Infers trees for a small built-in dataset in the default and extended "
       "modes and with support, at one and several threads, and compares them "
       "with the trees this release is known to produce."},
      {"completion", 0, "qraft completion bash|zsh|fish",
       "print a shell completion script",
       "Prints a completion script. For bash: qraft completion bash > "
       "~/.local/share/bash-completion/completions/qraft"},
      {"help", 0, "qraft help [COMMAND|TOPIC]", "show help on a command or topic",
       "Topics: options, modes, guide, support, ablations, test, output."},
  };
  return table;
}

const CommandSpec* find_command(const std::string& name) {
  for (const CommandSpec& c : command_table()) {
    if (name == c.name) return &c;
  }
  return nullptr;
}

std::string suggest(const std::string& word,
                    const std::vector<std::string>& candidates) {
  const auto distance = [](const std::string& a, const std::string& b) {
    std::vector<std::size_t> row(b.size() + 1);
    for (std::size_t j = 0; j <= b.size(); ++j) row[j] = j;
    for (std::size_t i = 1; i <= a.size(); ++i) {
      std::size_t diagonal = row[0];
      row[0] = i;
      for (std::size_t j = 1; j <= b.size(); ++j) {
        const std::size_t up = row[j];
        row[j] = std::min({row[j] + 1, row[j - 1] + 1,
                           diagonal + (a[i - 1] == b[j - 1] ? 0u : 1u)});
        diagonal = up;
      }
    }
    return row[b.size()];
  };
  std::string best;
  std::size_t best_distance = word.size() / 3 + 2;
  for (const std::string& c : candidates) {
    const std::size_t d = distance(word, c);
    if (d < best_distance) {
      best_distance = d;
      best = c;
    }
  }
  return best;
}

void parse_options(int begin, int argc, char** argv, unsigned command,
                   CliState& state) {
  const std::vector<OptionSpec>& table = option_table();
  bool options_done = false;
  for (int i = begin; i < argc; ++i) {
    std::string arg = argv[i];
    if (options_done || arg.size() < 2 || arg[0] != '-') {
      state.positionals.push_back(arg);
      continue;
    }
    if (arg == "--") {
      options_done = true;
      continue;
    }
    std::string inline_value;
    bool has_inline = false;
    if (arg.rfind("--", 0) == 0) {
      const std::size_t eq = arg.find('=');
      if (eq != std::string::npos) {
        inline_value = arg.substr(eq + 1);
        arg = arg.substr(0, eq);
        has_inline = true;
      }
    }
    const OptionSpec* spec = nullptr;
    for (const OptionSpec& o : table) {
      if (arg == o.long_name || (o.short_name != nullptr && arg == o.short_name)) {
        spec = &o;
        break;
      }
    }
    if (spec == nullptr) {
      std::vector<std::string> names;
      for (const OptionSpec& o : table) {
        if ((o.commands & command) != 0) names.push_back(o.long_name);
      }
      const std::string close = suggest(arg, names);
      throw UsageError("unknown argument: " + arg +
                       (close.empty() ? "" : " (did you mean " + close + "?)"));
    }
    if ((spec->commands & command) == 0) {
      const CommandSpec* owner = nullptr;
      for (const CommandSpec& c : command_table()) {
        if ((c.mask & command) != 0) owner = &c;
      }
      throw UsageError(std::string(spec->long_name) +
                       " is not an option of 'qraft " +
                       (owner != nullptr ? owner->name : "") + "'");
    }
    std::string value;
    if (spec->metavar != nullptr) {
      if (has_inline) {
        value = inline_value;
      } else {
        if (i + 1 >= argc) throw UsageError("missing value for " + arg);
        value = argv[++i];
      }
      if (spec->choices != nullptr) {
        std::vector<std::string> allowed;
        std::stringstream in(spec->choices);
        for (std::string item; std::getline(in, item, '|');) {
          allowed.push_back(item);
        }
        if (std::find(allowed.begin(), allowed.end(), value) == allowed.end()) {
          throw UsageError("invalid value for " + std::string(spec->long_name) +
                           ": " + value + " (expected " + spec->choices + ")");
        }
      }
    } else if (has_inline) {
      throw UsageError(std::string(spec->long_name) + " takes no value");
    }
    spec->apply(state, value);
  }
}

namespace {

// Word-wraps `text` to `width` columns.
std::vector<std::string> wrap(const std::string& text, std::size_t width) {
  std::vector<std::string> lines;
  std::stringstream in(text);
  std::string word;
  std::string line;
  while (in >> word) {
    if (!line.empty() && line.size() + 1 + word.size() > width) {
      lines.push_back(line);
      line.clear();
    }
    line += (line.empty() ? "" : " ") + word;
  }
  if (!line.empty()) lines.push_back(line);
  return lines;
}

std::string option_left(const OptionSpec& o) {
  std::string left = o.short_name != nullptr
                         ? std::string("  ") + o.short_name + ", "
                         : std::string("      ");
  left += o.long_name;
  if (o.metavar != nullptr) left += std::string(" ") + o.metavar;
  return left;
}

}  // namespace

std::string render_options(unsigned command, bool common_only,
                           const char* group) {
  constexpr std::size_t kColumn = 28;
  constexpr std::size_t kWidth = 80;
  std::string out;
  for (const char* g : kGroups) {
    if (group != nullptr && std::string(g) != group) continue;
    std::string section;
    for (const OptionSpec& o : option_table()) {
      if (std::string(o.group) != g || (o.commands & command) == 0) continue;
      if (common_only && !o.common) continue;
      const std::string left = option_left(o);
      const std::vector<std::string> help = wrap(o.help, kWidth - kColumn);
      if (left.size() + 2 <= kColumn) {
        section += left + std::string(kColumn - left.size(), ' ') +
                   (help.empty() ? "" : help[0]) + "\n";
      } else {
        section += left + "\n" + std::string(kColumn, ' ') +
                   (help.empty() ? "" : help[0]) + "\n";
      }
      for (std::size_t k = 1; k < help.size(); ++k) {
        section += std::string(kColumn, ' ') + help[k] + "\n";
      }
    }
    if (section.empty()) continue;
    const char* heading = group_heading(g);
    out += "\n" + std::string(heading != nullptr ? heading : g) +
           (heading != nullptr ? "" : ":") + "\n" + section;
  }
  return out;
}

std::string completion_script(const std::string& shell) {
  std::vector<std::string> commands;
  for (const CommandSpec& c : command_table()) commands.push_back(c.name);
  const auto join = [](const std::vector<std::string>& items) {
    std::string out;
    for (const std::string& s : items) out += (out.empty() ? "" : " ") + s;
    return out;
  };
  const auto options_of = [](unsigned mask) {
    std::vector<std::string> names;
    for (const OptionSpec& o : option_table()) {
      if ((o.commands & mask) == 0) continue;
      names.push_back(o.long_name);
      if (o.short_name != nullptr) names.push_back(o.short_name);
    }
    return names;
  };
  if (shell == "bash" || shell == "zsh") {
    std::string s;
    if (shell == "zsh") s += "autoload -U +X bashcompinit && bashcompinit\n";
    s += "# qraft completion for " + shell + "\n_qraft() {\n";
    s += "  local cur=\"${COMP_WORDS[COMP_CWORD]}\" "
         "prev=\"${COMP_WORDS[COMP_CWORD-1]}\" cmd=infer\n";
    s += "  case \" " + join(commands) + " \" in *\" ${COMP_WORDS[1]} \"*) "
         "cmd=\"${COMP_WORDS[1]}\" ;; esac\n";
    s += "  case \"$prev\" in\n";
    for (const OptionSpec& o : option_table()) {
      if (o.metavar == nullptr) continue;
      std::string names = o.long_name;
      if (o.short_name != nullptr) names += std::string("|") + o.short_name;
      if (o.choices != nullptr) {
        std::string words = o.choices;
        std::replace(words.begin(), words.end(), '|', ' ');
        s += "    " + names + ") COMPREPLY=($(compgen -W \"" + words +
             "\" -- \"$cur\")); return ;;\n";
      } else if (std::string(o.metavar) == "PATH") {
        s += "    " + names + ") COMPREPLY=($(compgen -f -- \"$cur\")); return ;;\n";
      } else {
        s += "    " + names + ") return ;;\n";
      }
    }
    s += "    completion) COMPREPLY=($(compgen -W \"bash zsh fish\" -- \"$cur\")); return ;;\n";
    s += "    help) COMPREPLY=($(compgen -W \"" + join(commands) +
         " options modes guide support ablations test output\" -- \"$cur\")); return ;;\n";
    s += "  esac\n";
    s += "  if [[ $COMP_CWORD -eq 1 && $cur != -* ]]; then\n";
    s += "    COMPREPLY=($(compgen -W \"" + join(commands) + "\" -- \"$cur\") "
         "$(compgen -f -- \"$cur\")); return\n  fi\n";
    s += "  case \"$cmd\" in\n";
    for (const CommandSpec& c : command_table()) {
      if (c.mask == 0) continue;
      s += "    " + std::string(c.name) + ") COMPREPLY=($(compgen -W \"" +
           join(options_of(c.mask)) + "\" -- \"$cur\")) ;;\n";
    }
    s += "  esac\n";
    s += "  [[ $cur != -* ]] && COMPREPLY+=($(compgen -f -- \"$cur\"))\n";
    s += "}\ncomplete -F _qraft qraft\n";
    return s;
  }
  if (shell == "fish") {
    std::string s = "# qraft completion for fish\n";
    s += "complete -c qraft -f\n";
    for (const CommandSpec& c : command_table()) {
      s += "complete -c qraft -n __fish_use_subcommand -a " + std::string(c.name) +
           " -d '" + c.summary + "'\n";
    }
    for (const OptionSpec& o : option_table()) {
      std::string line = "complete -c qraft -l " + std::string(o.long_name + 2);
      if (o.short_name != nullptr) line += std::string(" -s ") + (o.short_name + 1);
      if (o.metavar != nullptr) {
        line += " -r";
        if (o.choices != nullptr) {
          std::string words = o.choices;
          std::replace(words.begin(), words.end(), '|', ' ');
          line += " -a '" + words + "'";
        } else if (std::string(o.metavar) == "PATH") {
          line += " -F";
        }
      }
      std::string desc = o.help;
      if (desc.size() > 60) desc = desc.substr(0, 57) + "...";
      std::replace(desc.begin(), desc.end(), '\'', ' ');
      s += line + " -d '" + desc + "'\n";
    }
    return s;
  }
  throw UsageError("unknown shell: " + shell + " (expected bash, zsh or fish)");
}

std::string man_page(const std::string& version) {
  const auto roff = [](std::string text) {
    std::string out;
    for (char ch : text) {
      if (ch == '\\') out += "\\e";
      else if (ch == '-') out += "\\-";
      else out += ch;
    }
    if (!out.empty() && (out[0] == '.' || out[0] == '\'')) out = "\\&" + out;
    return out;
  };
  std::string s;
  s += ".TH QRAFT 1 \"\" \"qraft " + version + "\" \"User Commands\"\n";
  s += ".SH NAME\nqraft \\- species trees from gene trees by quartet regrafting\n";
  s += ".SH SYNOPSIS\n";
  for (const CommandSpec& c : command_table()) {
    s += ".B " + roff(c.usage) + "\n.br\n";
  }
  s += ".SH DESCRIPTION\n" +
       roff("QRAFT infers a species tree from gene trees by maximizing quartet "
            "support. It starts from the majority-rule extended consensus of "
            "the gene trees, divides the taxa recursively along well-supported "
            "bipartitions, joins the solved sides back by contextual "
            "regrafting and refines the result by nearest-neighbour "
            "interchanges, with every quartet count exact. The tree does not "
            "depend on the number of threads.") +
       "\n";
  s += ".SH COMMANDS\n";
  for (const CommandSpec& c : command_table()) {
    s += ".TP\n.B " + std::string(c.name) + "\n" + roff(c.description) + "\n";
  }
  for (const char* g : kGroups) {
    bool any = false;
    for (const OptionSpec& o : option_table()) {
      if (std::string(o.group) != g) continue;
      if (!any) {
        s += ".SH " + roff(std::string(g)) + "\n";
        any = true;
      }
      std::string head = o.short_name != nullptr
                             ? "\\fB" + roff(o.short_name) + "\\fR, "
                             : std::string();
      head += "\\fB" + roff(o.long_name) + "\\fR";
      if (o.metavar != nullptr) head += " \\fI" + roff(o.metavar) + "\\fR";
      s += ".TP\n" + head + "\n" + roff(o.help) + "\n";
    }
  }
  s += ".SH EXAMPLES\n.nf\n";
  s += roff("qraft -i genes.tre -o species.tre\n"
            "qraft -i genes.tre -o species.tre -t 8 --mode extended\n"
            "qraft -i genes.tre -o species.tre --guide-tree backbone.tre\n"
            "qraft -i genes.tre -o species.tre --support localpp\n"
            "qraft support -i genes.tre --tree species.tre -o annotated.tre\n"
            "qraft check genes.tre\n");
  s += ".fi\n";
  s += ".SH EXIT STATUS\n" +
       roff("0 on success; 1 when an input cannot be read or used, or the run "
            "fails; 2 for a command-line error; 3 when memory runs out.") +
       "\n";
  s += ".SH ENVIRONMENT\n.TP\n.B NO_COLOR\n" +
       roff("when set, messages are never coloured.") + "\n";
  s += ".SH LICENSE\n" +
       roff("GNU Affero General Public License, version 3 or later.") + "\n";
  return s;
}

}  // namespace qraft_cli
