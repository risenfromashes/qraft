#include "internal.hpp"

#include <memory>
#include <thread>

namespace qraft::internal {

namespace {
// The executor shared by the whole process. The inference run gate replaces it
// only while no other run uses it.
std::unique_ptr<tf::Executor> g_executor;
}

/// Replaces the shared executor with one of options.threads workers (every
/// hardware thread when 0) and tunes the allocator for that count.
int configure_threads(const Options& options) {
  int nw = options.threads > 0 ? options.threads
                               : static_cast<int>(std::thread::hardware_concurrency());
  if (nw < 1) nw = 1;
  configure_allocator(nw);
  g_executor = std::make_unique<tf::Executor>(nw);
  return nw;
}

/// Returns the shared executor, creating one with a worker per hardware
/// thread on first use.
tf::Executor& get_executor() {
  if (!g_executor) {
    g_executor = std::make_unique<tf::Executor>();
  }
  return *g_executor;
}

int get_num_workers() {
  return static_cast<int>(get_executor().num_workers());
}

/// Runs `taskflow` to completion. A worker thread co-runs it, executing other
/// tasks while it waits, so nested flows cannot deadlock the pool.
void run_inner(tf::Taskflow& taskflow) {
  auto& exec = get_executor();
  if (exec.this_worker_id() >= 0) {
    exec.corun(taskflow);
  } else {
    exec.run(taskflow).wait();
  }
}

}  // namespace qraft::internal
