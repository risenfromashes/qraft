#include "gene_input.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <new>
#include <stdexcept>
#include <string_view>
#include <system_error>
#include <thread>

#if defined(_WIN32)
#include <io.h>
#define QRAFT_ACCESS _access
#define QRAFT_R_OK 4
#define QRAFT_W_OK 2
#define QRAFT_WX_OK 2
#else
#include <unistd.h>
#define QRAFT_ACCESS access
#define QRAFT_R_OK R_OK
#define QRAFT_W_OK W_OK
#define QRAFT_WX_OK (W_OK | X_OK)
#endif

namespace qraft_cli {
namespace {

bool blank(std::string_view line) {
  return std::all_of(line.begin(), line.end(), [](char ch) {
    return std::isspace(static_cast<unsigned char>(ch)) != 0;
  });
}

// The genes file is read and parsed in blocks of this many bytes, so the text
// is never held in full next to the parsed trees.
constexpr std::size_t kGeneBlockBytes = std::size_t{256} << 20;
// A block is read by one thread per this many bytes, up to `threads`.
constexpr std::size_t kGeneReadBytesPerThread = std::size_t{8} << 20;

// Runs work(t) for t in [0, count) on `count` threads, the calling thread
// being thread 0.
template <class F>
void run_on_threads(int count, F&& work) {
  if (count <= 1) {
    work(0);
    return;
  }
  // The pool joins its threads when it goes out of scope, on unwinding too,
  // so a failed launch cannot outlive the work. (libc++, the C++ library of
  // macOS, has no std::jthread by default.)
  struct Pool {
    std::vector<std::thread> threads;
    ~Pool() {
      for (std::thread& thread : threads) thread.join();
    }
  } pool;
  pool.threads.reserve(static_cast<std::size_t>(count - 1));
  for (int t = 1; t < count; ++t) {
    pool.threads.emplace_back([&work, t] { work(t); });
  }
  work(0);
}

// Reads bytes [offset, offset + size) of `path` into `out` on up to
// `threads` threads, each with its own stream.
void read_file_range(const std::string& path,
                     std::uint64_t offset,
                     char* out,
                     std::size_t size,
                     int threads) {
  const int readers = static_cast<int>(std::clamp<std::size_t>(
      size / kGeneReadBytesPerThread, 1, static_cast<std::size_t>(threads)));
  std::vector<char> ok(static_cast<std::size_t>(readers), 0);
  run_on_threads(readers, [&](int t) {
    const std::size_t first = size * static_cast<std::size_t>(t) /
                              static_cast<std::size_t>(readers);
    const std::size_t last = size * static_cast<std::size_t>(t + 1) /
                             static_cast<std::size_t>(readers);
    std::ifstream in(path, std::ios::binary);
    if (!in) return;
    in.seekg(static_cast<std::streamoff>(offset + first));
    in.read(out + first, static_cast<std::streamsize>(last - first));
    ok[static_cast<std::size_t>(t)] =
        in.gcount() == static_cast<std::streamsize>(last - first) ? 1 : 0;
  });
  for (char success : ok) {
    if (success == 0) {
      throw std::runtime_error("cannot read the gene trees file '" + path + "'");
    }
  }
}

// The non-blank lines of text[0, size), which ends at a line end, as views
// in order; found on up to `threads` threads.
std::vector<std::string_view> split_gene_lines(const char* text,
                                               std::size_t size,
                                               int threads) {
  const int parts = static_cast<int>(std::clamp<std::size_t>(
      size / kGeneReadBytesPerThread, 1, static_cast<std::size_t>(threads)));
  // Part t holds the lines that start in its byte range.
  std::vector<std::vector<std::string_view>> found(
      static_cast<std::size_t>(parts));
  run_on_threads(parts, [&](int t) {
    std::size_t pos = size * static_cast<std::size_t>(t) /
                      static_cast<std::size_t>(parts);
    const std::size_t end = size * static_cast<std::size_t>(t + 1) /
                            static_cast<std::size_t>(parts);
    // A line starts at 0 or after a line end.
    if (t > 0) {
      while (pos < end && text[pos - 1] != '\n') ++pos;
    }
    while (pos < end) {
      const void* nl = std::memchr(text + pos, '\n', size - pos);
      const std::size_t stop =
          nl == nullptr ? size
                        : static_cast<std::size_t>(
                              static_cast<const char*>(nl) - text);
      const std::string_view line(text + pos, stop - pos);
      if (!blank(line)) found[static_cast<std::size_t>(t)].push_back(line);
      pos = stop + 1;
    }
  });
  std::size_t count = 0;
  for (const auto& part : found) count += part.size();
  std::vector<std::string_view> lines;
  lines.reserve(count);
  for (const auto& part : found) lines.insert(lines.end(), part.begin(), part.end());
  return lines;
}

// A line that failed to parse: its index in the block's lines and the
// parser's message.
struct LineFailure {
  std::size_t index = 0;
  std::string message;
};

// Parses `lines` into `genes`, appended in order, on up to `threads`
// threads. Returns the first failing line, whichever thread saw it, so the
// caller can name its line number.
bool parse_gene_lines(const std::vector<std::string_view>& lines,
                      int threads,
                      std::vector<qraft::GeneTree>& genes,
                      LineFailure& failure) {
  const std::size_t base = genes.size();
  genes.resize(base + lines.size());
  std::vector<std::exception_ptr> failures(lines.size());
  std::atomic<std::size_t> next{0};
  const int workers = static_cast<int>(std::clamp<std::size_t>(
      lines.size(), 1, static_cast<std::size_t>(threads)));
  run_on_threads(workers, [&](int) {
    for (std::size_t i; (i = next.fetch_add(1, std::memory_order_relaxed)) <
                        lines.size();) {
      try {
        genes[base + i] = qraft::parse_newick(lines[i]);
      } catch (...) {
        failures[i] = std::current_exception();
      }
    }
  });
  for (std::size_t i = 0; i < failures.size(); ++i) {
    if (failures[i] == nullptr) continue;
    failure.index = i;
    try {
      std::rethrow_exception(failures[i]);
    } catch (const std::bad_alloc&) {
      throw;
    } catch (const std::exception& ex) {
      failure.message = ex.what();
    } catch (...) {
      failure.message = "unknown parse error";
    }
    return false;
  }
  return true;
}

std::string error_text(int code) {
  return std::system_category().message(code);
}

}  // namespace

void require_readable(const std::string& path, const std::string& what,
                      bool nonempty) {
  std::error_code ec;
  const auto status = std::filesystem::status(path, ec);
  if (ec || !std::filesystem::exists(status)) {
    throw std::runtime_error("cannot read the " + what + " '" + path +
                             "': no such file");
  }
  if (std::filesystem::is_directory(status)) {
    throw std::runtime_error("cannot read the " + what + " '" + path +
                             "': it is a directory");
  }
  if (!std::filesystem::is_regular_file(status)) {
    // A pipe or another special file (such as <(...)) is read once, by the
    // run; opening it here could block or take its data.
    if (QRAFT_ACCESS(path.c_str(), QRAFT_R_OK) != 0) {
      throw std::runtime_error("cannot read the " + what + " '" + path +
                               "': " + error_text(errno));
    }
    return;
  }
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("cannot read the " + what + " '" + path +
                             "': " + error_text(errno));
  }
  if (nonempty && std::filesystem::file_size(path, ec) == 0 && !ec) {
    throw std::runtime_error("the " + what + " '" + path + "' is empty");
  }
}

void require_regular_gene_file(const std::string& path) {
  std::error_code ec;
  if (!std::filesystem::is_regular_file(std::filesystem::status(path, ec)) ||
      ec) {
    throw std::runtime_error(
        "the gene trees file '" + path +
        "' is a pipe or special file; it is read in parallel blocks, so "
        "save it to a regular file first");
  }
}

void check_tree_file(const std::string& path, const std::string& what) {
  std::error_code ec;
  if (!std::filesystem::is_regular_file(std::filesystem::status(path, ec)) ||
      ec) {
    return;
  }
  std::ifstream in(path, std::ios::binary);
  const std::string text((std::istreambuf_iterator<char>(in)),
                         std::istreambuf_iterator<char>());
  if (in.bad()) {
    throw std::runtime_error("cannot read the " + what + " '" + path + "'");
  }
  try {
    (void)qraft::parse_newick(text);
  } catch (const std::bad_alloc&) {
    throw;
  } catch (const std::exception& ex) {
    throw std::runtime_error("the " + what + " '" + path +
                             "' is not a valid Newick tree: " + ex.what());
  }
}

void require_writable(const std::string& path, const std::string& what) {
  namespace fs = std::filesystem;
  std::error_code ec;
  fs::path target(path);
  // Writing through a symbolic link to a missing file creates its target,
  // so the target's folder is the one that must exist.
  for (int hops = 0; hops < 40 && fs::is_symlink(fs::symlink_status(target, ec));
       ++hops) {
    const fs::path link = fs::read_symlink(target, ec);
    if (ec) break;
    target = link.is_absolute() ? link : target.parent_path() / link;
  }
  const auto status = fs::status(target, ec);
  if (!ec && fs::exists(status)) {
    if (fs::is_directory(status)) {
      throw std::runtime_error("cannot write the " + what + " '" + path +
                               "': it is a directory");
    }
    if (QRAFT_ACCESS(path.c_str(), QRAFT_W_OK) != 0) {
      throw std::runtime_error("cannot write the " + what + " '" + path +
                               "': " + error_text(errno));
    }
    return;
  }
  fs::path dir = target.parent_path();
  if (dir.empty()) dir = ".";
  const auto dir_status = fs::status(dir, ec);
  if (ec || !fs::is_directory(dir_status)) {
    throw std::runtime_error("cannot write the " + what + " '" + path +
                             "': the folder '" + dir.string() +
                             "' does not exist");
  }
  if (QRAFT_ACCESS(dir.string().c_str(), QRAFT_WX_OK) != 0) {
    throw std::runtime_error("cannot write the " + what + " '" + path +
                             "': the folder '" + dir.string() +
                             "' is not writable");
  }
}

std::vector<qraft::GeneTree> read_gene_trees(const std::string& path,
                                             int threads) {
  std::ifstream probe(path, std::ios::binary | std::ios::ate);
  if (!probe) {
    throw std::runtime_error("cannot read the gene trees file '" + path + "'");
  }
  const std::streamoff file_size = probe.tellg();
  probe.close();
  if (file_size < 0) {
    throw std::runtime_error("cannot read the gene trees file '" + path + "'");
  }
  threads = std::max(1, threads);

  std::vector<qraft::GeneTree> genes;
  std::string block;
  std::size_t carried = 0;
  std::uint64_t offset = 0;
  // Line ends before the current block, for error messages.
  std::uint64_t lines_before = 0;
  const std::uint64_t total = static_cast<std::uint64_t>(file_size);
  while (true) {
    const std::size_t fresh = static_cast<std::size_t>(std::min<std::uint64_t>(
        kGeneBlockBytes, total - offset));
    // The carried partial line is at the front of the block.
    block.resize(carried + fresh);
    if (fresh > 0) {
      read_file_range(path, offset, block.data() + carried, fresh, threads);
    }
    offset += fresh;
    // Lines up to the last line end are complete; at the end of the file the
    // rest is the last line.
    const bool last = offset >= total;
    std::size_t complete = block.size();
    if (!last) {
      const std::size_t nl = block.rfind('\n');
      complete = nl == std::string::npos ? 0 : nl + 1;
    }
    if (complete > 0) {
      const std::vector<std::string_view> lines =
          split_gene_lines(block.data(), complete, threads);
      LineFailure failure;
      if (!parse_gene_lines(lines, threads, genes, failure)) {
        const char* start = lines[failure.index].data();
        const std::uint64_t line =
            lines_before + 1 +
            static_cast<std::uint64_t>(std::count(
                static_cast<const char*>(block.data()), start, '\n'));
        throw std::runtime_error("the gene trees file '" + path + "', line " +
                                 std::to_string(line) +
                                 ", is not a valid Newick tree: " +
                                 failure.message);
      }
      lines_before += static_cast<std::uint64_t>(
          std::count(block.data(), block.data() + complete, '\n'));
    }
    if (last) break;
    carried = block.size() - complete;
    std::memmove(block.data(), block.data() + complete, carried);
  }
  if (genes.empty()) {
    throw std::runtime_error("the gene trees file '" + path +
                             "' holds no Newick trees (only blank lines)");
  }
  if (genes.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error("gene count exceeds the integer-indexed domain");
  }
  return genes;
}

}  // namespace qraft_cli
