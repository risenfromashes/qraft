function(qraft_prepare_deterministic_taskflow source_root destination_root)
  set(source_executor "${source_root}/taskflow/core/executor.hpp")
  # Hash and patch the header with LF line endings, so a checkout with CRLF
  # line endings (the Git for Windows default) gives the same overlay.
  file(READ "${source_executor}" source_contents)
  string(REPLACE "\r\n" "\n" source_contents "${source_contents}")
  string(SHA256 source_executor_sha256 "${source_contents}")
  if(NOT source_executor_sha256 STREQUAL
     "eeb9ba321e30aeb24615b4e24344aa83b77804a8bf40a340a6df7f6bf6b507d6")
    message(FATAL_ERROR
      "Taskflow executor changed; re-audit the deterministic worker patch")
  endif()

  set(patched_executor "${destination_root}/taskflow/core/executor.hpp")
  set(expected_patched_sha256
    "6bdc08848ef9a1cbb7b8c4fa5404aa8ab32e58275b4189c33a54eb1948e68928")
  set(materialize TRUE)
  if(EXISTS "${patched_executor}")
    file(SHA256 "${patched_executor}" existing_patched_sha256)
    if(existing_patched_sha256 STREQUAL expected_patched_sha256)
      set(materialize FALSE)
    endif()
  endif()

  if(materialize)
  file(REMOVE_RECURSE "${destination_root}/taskflow")
  file(MAKE_DIRECTORY "${destination_root}")
  file(COPY "${source_root}/taskflow" DESTINATION "${destination_root}")
  set(contents "${source_contents}")

  string(REPLACE
    "#pragma once\n"
    "#pragma once\n\n#define QRAFT_TASKFLOW_DETERMINISTIC_WORKER_LOOKUP 1\n"
    contents "${contents}")

  string(REPLACE
    "  std::unordered_map<std::thread::id, Worker*> _t2w;"
    [=[  // Production uses one process-global executor. The executor pointer
  // also makes independent test executors safe on their own worker threads.
  inline static thread_local Executor* _tls_executor {nullptr};
  inline static thread_local Worker* _tls_worker {nullptr};]=]
    contents "${contents}")

  string(REPLACE
    [=[inline Worker* Executor::this_worker() {
  auto itr = _t2w.find(std::this_thread::get_id());
  return itr == _t2w.end() ? nullptr : itr->second;
}

// Function: this_worker_id
inline int Executor::this_worker_id() const {
  auto i = _t2w.find(std::this_thread::get_id());
  return i == _t2w.end() ? -1 : static_cast<int>(i->second->_id);
}]=]
    [=[inline Worker* Executor::this_worker() {
  return _tls_executor == this ? _tls_worker : nullptr;
}

// Function: this_worker_id
inline int Executor::this_worker_id() const {
  return _tls_executor == this && _tls_worker != nullptr
           ? static_cast<int>(_tls_worker->_id)
           : -1;
}]=]
    contents "${contents}")

  string(CONCAT source_worker_init
    "      auto& worker = _workers[id];\n"
    "  \n"
    "      worker._id = id;")
  string(CONCAT deterministic_worker_init
    "      auto& worker = _workers[id];\n"
    "      _tls_executor = this;\n"
    "      _tls_worker = &worker;\n"
    "  \n"
    "      worker._id = id;")
  string(REPLACE "${source_worker_init}" "${deterministic_worker_init}"
    contents "${contents}")

  string(REPLACE
    [=[      if(wif) {
        wif->scheduler_epilogue(worker, ptr);
      }

    });]=]
    [=[      if(wif) {
        wif->scheduler_epilogue(worker, ptr);
      }

      _tls_worker = nullptr;
      _tls_executor = nullptr;

    });]=]
    contents "${contents}")

  string(CONCAT source_worker_registration
    "    \n"
    "    // We avoid using thread-local storage to track the mapping between a thread\n"
    "    // and its corresponding worker in an executor. On Windows, thread-local\n"
    "    // storage can be unreliable in certain situations (see issue #727).\n"
    "    //\n"
    "    // Instead, we maintain a per-executor mapping from threads to workers.\n"
    "    // This approach has an additional advantage: according to the C++ Standard,\n"
    "    // std::thread::id uniquely identifies a thread object. Therefore, once the map\n"
    "    // returns a valid worker, we can be certain that the worker belongs to this\n"
    "    // executor. This eliminates the need for additional executor validation \n"
    "    // required by using thread-local storage.\n"
    "    //\n"
    "    // Example:\n"
    "    //\n"
    "    //   Worker* w = this_worker();\n"
    "    //   // Using thread-local storage, we would need additional executor validation:\n"
    "    //   if (w == nullptr || w->_executor != this) { /* caller is not a worker of this executor */ }\n"
    "    //\n"
    "    //   // Using per-executor mapping, it suffices to check:\n"
    "    //   if (w == nullptr) { /* caller is not a worker of this executor */ }\n"
    "    //\n"
    "    _t2w.emplace(_workers[id]._thread.get_id(), &_workers[id]);")
  string(REPLACE "${source_worker_registration}" "" contents "${contents}")

  if(contents MATCHES "_t2w" OR
     NOT contents MATCHES "QRAFT_TASKFLOW_DETERMINISTIC_WORKER_LOOKUP" OR
     NOT contents MATCHES "thread_local Executor\\* _tls_executor" OR
     NOT contents MATCHES "thread_local Worker\\* _tls_worker")
    message(FATAL_ERROR
      "failed to materialize deterministic Taskflow worker lookup")
  endif()

  file(WRITE "${patched_executor}" "${contents}")
  endif()

  file(SHA256 "${patched_executor}" patched_sha256)
  if(NOT patched_sha256 STREQUAL expected_patched_sha256)
    message(FATAL_ERROR "deterministic Taskflow overlay hash mismatch")
  endif()
  set(QRAFT_TASKFLOW_INCLUDE_ROOT "${destination_root}" PARENT_SCOPE)
  set(QRAFT_TASKFLOW_EXECUTOR "${patched_executor}" PARENT_SCOPE)
endfunction()
