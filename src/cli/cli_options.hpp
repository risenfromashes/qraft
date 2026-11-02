#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include "qraft/options.hpp"

namespace qraft_cli {

/// Commands that take options, as a bit mask so an option can belong to
/// several. Adding a command means one entry in command_table() and its bit
/// on the options it takes.
enum Command : unsigned {
  kInfer = 1u << 0,
  kSupport = 1u << 1,
  kScore = 1u << 2,
  kCheck = 1u << 3,
};

enum class Verbosity { Quiet, Normal, Verbose };

/// Everything the command line sets.
struct CliState {
  std::string genes_path;
  std::string out_path;
  std::string tree_path;  ///< support and score: the species tree to read
  std::string json_path;  ///< check: the summary as JSON
  qraft::Options options;
  Verbosity verbosity = Verbosity::Normal;
  bool progress_given = false;  ///< --progress or --no-progress was given
  bool help = false;
  std::vector<std::string> positionals;
};

/// One command-line option. The same table drives the parser, the help
/// text, the man page and the shell completions.
struct OptionSpec {
  const char* long_name;   ///< "--genes"
  const char* short_name;  ///< "-i", or nullptr
  const char* metavar;     ///< "PATH" for an option with a value, else nullptr
  const char* choices;     ///< "a|b" when the value is one of these, else nullptr
  const char* group;       ///< help section
  unsigned commands;       ///< Command mask
  bool common;             ///< listed in the short help
  const char* help;
  void (*apply)(CliState&, const std::string&);
};

struct CommandSpec {
  const char* name;
  unsigned mask;  ///< its Command bit, or 0 for commands without options
  const char* usage;
  const char* summary;
  const char* description;
};

const std::vector<OptionSpec>& option_table();
const std::vector<CommandSpec>& command_table();
const CommandSpec* find_command(const std::string& name);

/// A command-line error, reported with a pointer to --help and exit status 2.
struct UsageError : std::invalid_argument {
  using std::invalid_argument::invalid_argument;
};

/// `value` as an int; throws UsageError naming `flag` unless all of it is a
/// whole number.
int parse_int(const std::string& flag, const std::string& value);

/// Applies argv[begin, argc) to `state` for `command`; throws UsageError.
void parse_options(int begin, int argc, char** argv, unsigned command,
                   CliState& state);

/// The closest of `candidates` to `word`, or "" when none is close.
std::string suggest(const std::string& word,
                    const std::vector<std::string>& candidates);

/// The options of `command` as help text, grouped; only the common ones when
/// `common_only`; only `group` when it is not null.
std::string render_options(unsigned command, bool common_only,
                           const char* group = nullptr);
std::string completion_script(const std::string& shell);
std::string man_page(const std::string& version);

}  // namespace qraft_cli
