#pragma once

#include <cstddef>
#include <exception>
#include <functional>
#include <stdexcept>

#if defined(__linux__)
#include <pthread.h>
#endif

namespace qraft::test {

#if defined(__linux__)
namespace detail {

struct SmallStackCall {
  const std::function<void()>* fn = nullptr;
  std::exception_ptr failure;
};

inline void* run_small_stack_call(void* argument) {
  auto* call = static_cast<SmallStackCall*>(argument);
  try {
    (*call->fn)();
  } catch (...) {
    call->failure = std::current_exception();
  }
  return nullptr;
}

}  // namespace detail
#endif

/// Runs `fn` on a new thread with a stack of `bytes` and rethrows its
/// exception. A recursion as deep as its input overflows such a stack and
/// crashes the test. Elsewhere than Linux, runs `fn` on the calling thread.
inline void run_with_stack(std::size_t bytes, const std::function<void()>& fn) {
#if defined(__linux__)
  detail::SmallStackCall call;
  call.fn = &fn;
  pthread_attr_t attributes;
  if (pthread_attr_init(&attributes) != 0) {
    throw std::runtime_error("failed to initialise thread attributes");
  }
  if (pthread_attr_setstacksize(&attributes, bytes) != 0) {
    pthread_attr_destroy(&attributes);
    throw std::runtime_error("failed to set a thread stack size");
  }
  pthread_t thread;
  const int created = pthread_create(&thread, &attributes,
                                     detail::run_small_stack_call, &call);
  pthread_attr_destroy(&attributes);
  if (created != 0) {
    throw std::runtime_error("failed to start a small-stack thread");
  }
  pthread_join(thread, nullptr);
  if (call.failure) std::rethrow_exception(call.failure);
#else
  (void)bytes;
  fn();
#endif
}

/// Gives threads created from now on, such as the workers of a new executor,
/// stacks of `bytes`. Returns false where this cannot be set.
inline bool set_default_thread_stack(std::size_t bytes) {
#if defined(__GLIBC__)
  pthread_attr_t attributes;
  if (pthread_getattr_default_np(&attributes) != 0) return false;
  const bool set = pthread_attr_setstacksize(&attributes, bytes) == 0 &&
                   pthread_setattr_default_np(&attributes) == 0;
  pthread_attr_destroy(&attributes);
  return set;
#else
  (void)bytes;
  return false;
#endif
}

}  // namespace qraft::test
