// Allocator hooks. Builds that link jemalloc, as the default static build
// always does, define QRAFT_LINKED_JEMALLOC; other builds use the C library's
// malloc, tuned here when it is glibc's. A jemalloc build does not reference
// glibc's malloc controls, since they would pull glibc's malloc into the
// static link.

#include "internal.hpp"

#include <algorithm>

#if defined(QRAFT_LINKED_JEMALLOC)
extern "C" int mallctl(const char* name, void* oldp, std::size_t* oldlenp,
                        void* newp, std::size_t newlen);
#elif defined(__GLIBC__)
#include <malloc.h>
#endif

namespace qraft::internal {

void configure_allocator(int workers) {
#if !defined(QRAFT_LINKED_JEMALLOC) && defined(__GLIBC__)
  // At most two arenas, and free memory above 4 MiB at the top of the heap
  // is returned to the system.
  mallopt(M_ARENA_MAX, std::max(1, std::min(workers, 2)));
  mallopt(M_TRIM_THRESHOLD, 4 * 1024 * 1024);
#else
  (void)workers;
#endif
}

void release_retained_allocator_pages() {
#if defined(QRAFT_LINKED_JEMALLOC)
  // 4096 is MALLCTL_ARENAS_ALL. A failed purge only leaves pages resident.
  (void)mallctl("arena.4096.purge", nullptr, nullptr, nullptr, 0);
#endif
}

void release_after_subproblem(int active_taxa) {
#if defined(QRAFT_LINKED_JEMALLOC)
  if (active_taxa >= 512) release_retained_allocator_pages();
#elif defined(__GLIBC__)
  if (active_taxa >= 128) malloc_trim(0);
#else
  (void)active_taxa;
#endif
}

}  // namespace qraft::internal
