#pragma once

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "internal.hpp"
#include "solver/solver.hpp"

namespace qraft::internal {

/// The stage a run is in, for the status file.
enum class RunStage : int {
  Startup = 0,
  StartingTree,
  Divide,
  Nni,
  Spr,
  Support,
  Done,
};

/// Notes the stage the run is in.
void checkpoint_note_stage(RunStage stage);

/// The outcome of one scored SPR candidate: whether its clade moved, the
/// edge it was cut at, and the two regraft edges it was joined at (in the
/// numbering of the two cut sides).
struct SprOutcome {
  int index = -1;
  bool moved = false;
  int edge_a = -1;
  int edge_b = -1;
  std::array<int, 4> joins{-1, -1, -1, -1};
};

/// The journaled outcomes of one SPR round, in candidate order. A complete
/// round was journaled before its nested NNI.
struct SprRoundLog {
  bool complete = false;
  std::vector<SprOutcome> outcomes;
};

/// Opens the run's journal and status file (both optional), owning the
/// process-wide checkpoint state for one run. The run must be exclusive.
/// `input_fingerprint` identifies the indexed input; a resumed run must match
/// it, the options and the program version.
class ScopedCheckpointRun {
 public:
  ScopedCheckpointRun(const Options& options, std::string input_fingerprint);
  ScopedCheckpointRun(const ScopedCheckpointRun&) = delete;
  ScopedCheckpointRun& operator=(const ScopedCheckpointRun&) = delete;
  /// Writes pending decisions and the final status; never throws.
  ~ScopedCheckpointRun();
  /// Marks the run finished: journals the end and writes the final status.
  void finish();

 private:
  bool active_ = false;
};

/// A hash of the gene trees as parsed: every node's label and neighbours, in
/// order. Computed only for journaling runs.
uint64_t checkpoint_genes_hash(const std::vector<GeneTree>& genes);

/// The canonical input fingerprint: gene count, taxon count, physical nodes,
/// the gene-tree hash, a hash of the taxon labels in id order, and hashes of
/// the starting-tree and guide-tree files.
std::string checkpoint_input_fingerprint(const Options& options,
                                         const TaxonRegistry& taxa,
                                         std::size_t gene_count,
                                         std::size_t physical_nodes,
                                         uint64_t genes_hash);

/// Whether the run journals decisions.
bool checkpoint_journaling();

/// The next invocation number of a pipeline step ('R' recursion, 'N' NNI,
/// 'S' SPR); invocations of one kind are numbered in call order. Returns -1
/// when the run does not journal.
int checkpoint_next_invocation(char kind);

/// The journaled left side of a subproblem's bipartition, as a bitset over
/// its `active` taxa, if any.
std::optional<Bitset> checkpoint_split(int invocation, uint64_t path,
                                       int active);
void checkpoint_record_split(int invocation, uint64_t path,
                             const Bitset& left_side);

/// The journaled regraft edges of a subproblem, if any.
std::optional<std::pair<RegraftEdge, RegraftEdge>> checkpoint_regraft(
    int invocation, uint64_t path);
void checkpoint_record_regraft(int invocation, uint64_t path,
                               const RegraftEdge& left,
                               const RegraftEdge& right);

/// The journaled NNI move groups of an invocation, in order: one move per
/// group, or one batch. Only a finished invocation has them: the groups of an
/// invocation the journal holds without its end are dropped when the journal
/// is read, and that invocation runs again from its start.
const std::vector<std::vector<NniMove>>* checkpoint_nni_groups(int invocation);
void checkpoint_record_nni_group(int invocation,
                                 const std::vector<NniMove>& group);
/// Whether the journal holds the end of NNI invocation `invocation`; and the
/// record of that end, written once the invocation has applied its last move.
bool checkpoint_nni_finished(int invocation);
void checkpoint_record_nni_end(int invocation);

/// The journaled outcomes of an SPR round, if any.
const SprRoundLog* checkpoint_spr_round(int invocation, int round);
void checkpoint_record_spr_outcome(int invocation, int round,
                                   const SprOutcome& outcome);
void checkpoint_record_spr_round_end(int invocation, int round);

/// Counts a decision taken from the journal instead of being computed.
void checkpoint_note_replayed();

}  // namespace qraft::internal
