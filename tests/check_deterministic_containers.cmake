# Fails when a library source names a hash container, whose iteration order
# could make construction or reduction order input-dependent, and checks that
# the Taskflow executor finds its worker through thread-local pointers rather
# than a thread-id hash map.

if(NOT DEFINED QRAFT_SOURCE_DIR)
  message(FATAL_ERROR "QRAFT_SOURCE_DIR is required")
endif()

file(GLOB_RECURSE production_files LIST_DIRECTORIES false
  "${QRAFT_SOURCE_DIR}/src/*.cpp"
  "${QRAFT_SOURCE_DIR}/src/*.hpp"
  "${QRAFT_SOURCE_DIR}/include/*.hpp")

foreach(path IN LISTS production_files)
  file(READ "${path}" contents)
  if(contents MATCHES "(unordered_map|unordered_set)")
    file(RELATIVE_PATH relative "${QRAFT_SOURCE_DIR}" "${path}")
    message(FATAL_ERROR
      "input-dependent hash-table type found in production source: ${relative}")
  endif()
endforeach()

if(NOT DEFINED QRAFT_TASKFLOW_EXECUTOR)
  message(FATAL_ERROR "QRAFT_TASKFLOW_EXECUTOR is required")
endif()
file(READ "${QRAFT_TASKFLOW_EXECUTOR}" taskflow_contents)
if(taskflow_contents MATCHES "_t2w" OR
   taskflow_contents MATCHES
   "unordered_map[ \\t\\r\\n]*<[^>]*thread::id[^>]*Worker")
  message(FATAL_ERROR
    "Taskflow worker lookup still depends on an unordered_map")
endif()
if(NOT taskflow_contents MATCHES
   "thread_local[ \\t\\r\\n]+Executor\\*[ \\t\\r\\n]+_tls_executor" OR
   NOT taskflow_contents MATCHES
   "thread_local[ \\t\\r\\n]+Worker\\*[ \\t\\r\\n]+_tls_worker")
  message(FATAL_ERROR
    "Taskflow deterministic TLS worker lookup contract is missing")
endif()

message(STATUS
  "production indexes and the active Taskflow worker lookup are deterministic")
