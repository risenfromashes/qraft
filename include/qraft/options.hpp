#pragma once

#include <string>

namespace qraft {

/// Progress reporting on standard error. The library reports only for On;
/// the command-line tool resolves Auto to On when standard error is a
/// terminal.
enum class ProgressMode { Off, Auto, On };

/// How the recursive divide ranks the bipartitions of the starting tree.
enum class DivideScore {
  Quartet,   ///< satisfied-minus-violated quartet count, W(00|11) - W(01|01)
  Balanced,  ///< most balanced edge, without quartet scoring
};

/// How contextual regrafting scores an attachment edge. Each internal vertex
/// of the side being attached to contributes the quartets with one taxon in
/// each of its three arms R0, R1, R2 and one on the other side.
enum class RegraftScore {
  Raw,         ///< unweighted quartet counts
  Normalized,  ///< each vertex's quartets divided by |R0||R1||R2||other side|
};

/// Branch support written on the internal edges of the output tree.
enum class SupportMode {
  None,     ///< topology only
  LocalPP,  ///< local posterior probability and length in coalescent units
};

/// The gene-tree quartet frequencies f1, f2, f3 of an edge's three
/// topologies, from the quartets with one taxon in each of its four arms.
/// The two agree when every gene holds every taxon and is binary.
enum class SupportFrequencies {
  Gene,     ///< each gene's three counts divided by their sum, then summed (ASTER; ASTRAL-III on binary genes)
  Quartet,  ///< counts summed over genes, divided by the quartets around the edge
};

/// The anchor-separation test applied to each proposed bipartition. A
/// proposal is accepted when its anchor separation exceeds the threshold t_k
/// for k genes; after `proposals` rejections the subproblem uses the trivial
/// bipartition (one taxon versus the rest). When genes miss taxa, the test
/// uses the ranks of the taxa along the anchor path in place of path lengths.
/// The test is skipped, accepting the top-ranked bipartition, when a gene
/// repeats a taxon.
struct BipartitionTestOptions {
  bool enabled = true;         ///< false accepts the top-ranked bipartition
  double coefficient = 0.05;   ///< C in t_k; finite and positive
  double crossover = 10000.0;  ///< K0 >= 0; t_k < 0 for k < K0, > 0 for k > K0
  int proposals = 1;           ///< bipartitions tested per subproblem, 1..16
  int anchor_pairs = 1;        ///< anchor pairs per bipartition, 1..4
};

/// Configuration of one inference.
struct Options {
  int threads = 0;  ///< worker threads; 0 uses every hardware thread
  ProgressMode progress = ProgressMode::Off;
  /// Newick file holding each taxon once, used in place of the MRE starting
  /// tree when set.
  std::string starting_tree_path;
  /// Newick guide tree on some or all taxa, when set. The starting tree
  /// loses its edges that conflict with the guide on the shared taxa and
  /// gains the guide bipartitions that keep its polytomies' arms whole.
  std::string guide_tree_path;
  /// Builds the MRE consensus without the clusters that conflict with the
  /// guide, when set; the guide then refines it as above. Needs
  /// guide_tree_path and no starting_tree_path.
  bool guide_seed = false;
  /// Writes the starting tree used, after the guide, when set.
  std::string write_starting_tree_path;
  std::string stats_json;  ///< writes the options and run statistics when set

  /// Recursive divide with contextual regrafting; without it the starting
  /// tree is refined directly.
  bool divide = true;
  DivideScore divide_score = DivideScore::Quartet;
  BipartitionTestOptions bipartition_test;
  RegraftScore regraft_score = RegraftScore::Raw;
  bool nni = true;  ///< gated NNI refinement of the final tree
  /// At most 20 NNI moves; false applies moves until none qualifies.
  bool nni_move_cap = true;
  /// The NNI support gate; false ignores the support of the replaced
  /// bipartition.
  bool nni_support_gate = true;
  /// Applies, per scoring pass, the set of NNI moves of largest total gain
  /// among sets whose gains add exactly, instead of one move at a time.
  bool nni_batch = false;
  /// SPR refinement after NNI, the second stage of extended mode: up to
  /// spr_rounds rounds, each ranking clades by their skeleton gain from a
  /// pass that reads at most spr_reads core vertices, then re-joining the
  /// first spr_candidates of them by the exact two-sided regraft and keeping
  /// every improving move, then NNI. A round without a move ends it.
  bool spr_refine = false;
  /// Core vertices read per pass, more only when more tie at the largest
  /// edge; positive.
  int spr_reads = 128;
  /// Fewest taxa of a ranked clade, the smaller side of its edge; at
  /// least 2.
  int spr_min_clade = 10;
  int spr_candidates = 150;  ///< clades verified per round; positive
  int spr_rounds = 2;        ///< positive

  /// Writes support:length on every internal edge whose ends are binary,
  /// from one more counting pass over the final tree; terminal edges and
  /// edges at a polytomy get none.
  SupportMode support = SupportMode::None;
  SupportFrequencies support_frequencies = SupportFrequencies::Gene;
  /// Writes the frequencies and the three posteriors of each edge as a
  /// bracketed label in place of the bare support.
  bool support_detail = false;

  /// Written as "qraft_version" in the statistics JSON when set (the command
  /// line sets it to the version and commit); it changes nothing else.
  std::string program_version;
};

}  // namespace qraft
