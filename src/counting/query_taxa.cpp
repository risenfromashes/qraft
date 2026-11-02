#include "counting/query_taxa.hpp"

#include <algorithm>

namespace qraft::internal {

int OwnerKeyMap::find(int key) const {
  if (key < 0) return -1;
  switch (mode) {
    case OwnerKeyMode::Missing:
      return -1;
    case OwnerKeyMode::Dense:
      return key < static_cast<int>(dense.size())
                 ? dense[static_cast<std::size_t>(key)]
                 : -1;
    case OwnerKeyMode::Sorted: {
      const auto it = std::lower_bound(
          sorted.begin(), sorted.end(), key,
          [](const std::pair<int, int>& entry, int value) {
            return entry.first < value;
          });
      return it != sorted.end() && it->first == key ? it->second : -1;
    }
  }
  return -1;
}

int OwnerKeyMap::inferred_owner_count() const {
  int count = 0;
  const auto note_owner = [&](int owner) {
    if (owner >= 0) count = std::max(count, owner + 1);
  };
  if (mode == OwnerKeyMode::Dense) {
    for (int owner : dense) note_owner(owner);
  } else if (mode == OwnerKeyMode::Sorted) {
    for (const auto& entry : sorted) note_owner(entry.second);
  }
  return count;
}

int owner_for_taxon(const OwnerProjection& owners, int taxon) {
  if (owners.taxon_keys.mode != OwnerKeyMode::Missing) {
    return owners.taxon_keys.find(taxon);
  }
  return taxon >= 0 && taxon < static_cast<int>(owners.taxon_to_owner.size())
             ? owners.taxon_to_owner[static_cast<std::size_t>(taxon)]
             : -1;
}

int owner_for_gene_leaf(const GeneIndex::Node& node,
                        const OwnerProjection& owners) {
  // A leaf with a local index maps only through local_keys; a miss there is
  // final, since the leaf carries no taxon to fall back on.
  const int local = node.local();
  const int owner = local >= 0 ? owners.local_keys.find(local)
                               : owner_for_taxon(owners, node.taxon());
  int owner_count = owners.owner_count;
  // An unset owner count is inferred from the maps.
  if (owner_count <= 0) {
    owner_count = std::max(owners.taxon_keys.inferred_owner_count(),
                           owners.local_keys.inferred_owner_count());
    for (int mapped : owners.taxon_to_owner) {
      owner_count = std::max(owner_count, mapped + 1);
    }
  }
  if (owner < 0 || owner >= owner_count) return -1;
  if (!owners.live_owner.empty() &&
      (owner >= static_cast<int>(owners.live_owner.size()) ||
       owners.live_owner[static_cast<std::size_t>(owner)] == 0)) {
    return -1;
  }
  return owner;
}

bool gene_has_owner_quartet(const GeneIndex& gene,
                            const OwnerProjection& owners) {
  int leaves = 0;
  for (const GeneIndex::Node& node : gene.nodes) {
    if (owner_for_gene_leaf(node, owners) >= 0 && ++leaves == 4) return true;
  }
  return false;
}

OwnerProjection make_owner_projection(const SolverContext& ctx,
                                      const TaxonIndex& local_by_taxon) {
  OwnerProjection owners;
  owners.owner_count = static_cast<int>(local_by_taxon.taxa.size());
  owners.live_owner.assign(static_cast<std::size_t>(owners.owner_count), 1);

  // A dense table is used only for a key domain of at most twice the owner
  // count, so the projection takes O(m) space for m owners.
  const auto linear_domain = [&](std::size_t domain) {
    return owners.owner_count > 0 &&
           domain <= static_cast<std::size_t>(owners.owner_count) * 2u;
  };
  std::size_t taxon_domain = context_taxa(ctx).size();
  if (!local_by_taxon.by_taxon.empty() &&
      local_by_taxon.by_taxon.back().first >= 0) {
    taxon_domain = std::max(
        taxon_domain,
        static_cast<std::size_t>(local_by_taxon.by_taxon.back().first) + 1u);
  }
  if (!local_by_taxon.by_taxon.empty()) {
    if (linear_domain(taxon_domain)) {
      owners.taxon_keys.mode = OwnerKeyMode::Dense;
      owners.taxon_keys.dense.assign(taxon_domain, -1);
      for (const auto& [taxon, owner] : local_by_taxon.by_taxon) {
        if (taxon >= 0 && taxon < static_cast<int>(taxon_domain)) {
          owners.taxon_keys.dense[static_cast<std::size_t>(taxon)] = owner;
        }
      }
    } else {
      owners.taxon_keys.mode = OwnerKeyMode::Sorted;
      owners.taxon_keys.sorted = local_by_taxon.by_taxon;
    }
  }

  if (ctx.local_taxa == nullptr) return owners;
  const std::vector<int>& local_taxa = *ctx.local_taxa;
  OwnerKeyMap& map = owners.local_keys;
  if (linear_domain(local_taxa.size())) {
    map.mode = OwnerKeyMode::Dense;
    map.dense.resize(local_taxa.size());
    for (std::size_t local = 0; local < local_taxa.size(); ++local) {
      map.dense[local] = owner_for_taxon(owners, local_taxa[local]);
    }
    return owners;
  }
  map.mode = OwnerKeyMode::Sorted;
  map.sorted.reserve(std::min(local_taxa.size(), local_by_taxon.taxa.size()));
  for (int local = 0; local < static_cast<int>(local_taxa.size()); ++local) {
    const int owner =
        owner_for_taxon(owners, local_taxa[static_cast<std::size_t>(local)]);
    if (owner >= 0) map.sorted.push_back({local, owner});
  }
  return owners;
}

}  // namespace qraft::internal
