#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "qraft/tree.hpp"

namespace qraft_cli {

/// Parses the non-blank lines of `path`, one Newick gene tree per line, in
/// input order on `threads` threads. A malformed line is reported with its
/// line number in the file and the parser's message.
std::vector<qraft::GeneTree> read_gene_trees(const std::string& path,
                                             int threads);

/// Throws std::runtime_error with a plain-language message unless `path`
/// exists, is not a directory and can be read; `what` names the file in the
/// message, such as "gene trees file". With `nonempty`, an empty regular file
/// is an error too. A pipe or other special file is not opened.
void require_readable(const std::string& path, const std::string& what,
                      bool nonempty);

/// Throws std::runtime_error unless `path` is a regular file: the gene trees
/// are read in parallel blocks, so a pipe cannot be used.
void require_regular_gene_file(const std::string& path);

/// Parses the whole of the regular file `path` as one Newick tree, as the run
/// reads a starting or guide tree, and throws std::runtime_error naming
/// `what` when it is malformed. A pipe or other special file is left to the
/// run.
void check_tree_file(const std::string& path, const std::string& what);

/// Throws std::runtime_error unless `path` can be written: an existing
/// writable file, or a new file in an existing folder that can be written and
/// searched, following symbolic links. Nothing is created.
void require_writable(const std::string& path, const std::string& what);

}  // namespace qraft_cli
