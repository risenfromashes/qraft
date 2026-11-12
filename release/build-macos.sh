#!/usr/bin/env bash
# Builds the macOS arm64 release archive on a Mac with Xcode's command line
# tools and CMake; the release workflow runs it on a GitHub macOS runner.
#
#   release/build-macos.sh [OUT_DIR]
#
# Apple clang, arm64, for macOS 12 or later. macOS has no static C library, so
# the binary links only the system's libSystem and libc++, dynamically, and
# uses the system allocator. The archive goes to OUT_DIR (default dist/); its
# times come from the commit date, as in release/build-release.sh.
set -euo pipefail
root=$(git -C "$(dirname "$0")/.." rev-parse --show-toplevel)
out=${1:-$root/dist}
mkdir -p "$out"
out=$(cd "$out" && pwd)
. "$root/release/names.sh"
. "$root/release/stage.sh"
target=macos-arm64
version=$(release_version "$root")
jobs=${JOBS:-$(sysctl -n hw.ncpu)}
SOURCE_DATE_EPOCH=$(git -C "$root" log -1 --format=%ct)
export SOURCE_DATE_EPOCH
git -C "$root" submodule update --init external/taskflow
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
build=$work/build
stage=$work/stage

cmake -S "$root" -B "$build" -DCMAKE_BUILD_TYPE=Release -DQRAFT_BUILD_TESTS=OFF \
  -DQRAFT_STATIC_EXECUTABLE=OFF -DQRAFT_USE_JEMALLOC=OFF \
  -DCMAKE_OSX_ARCHITECTURES=arm64 -DCMAKE_OSX_DEPLOYMENT_TARGET=12.0 \
  "-DCMAKE_CXX_FLAGS=-ffile-prefix-map=$root=. -ffile-prefix-map=$build=."
cmake --build "$build" --target qraft -j "$jobs"
qraft=$build/qraft
# Stripping invalidates the linker's ad-hoc signature, which arm64 macOS
# requires before it runs a binary; sign it again, ad hoc.
strip "$qraft"
codesign --force --sign - "$qraft"
codesign --verify --strict "$qraft"

# Only system libraries, and macOS 12 as the oldest system.
libraries=$(otool -L "$qraft" | tail -n +2 | awk '{print $1}')
for library in $libraries; do
  case $library in
    /usr/lib/libSystem.B.dylib | /usr/lib/libc++.1.dylib) ;;
    *) echo "build-macos.sh: unexpected library $library" >&2; exit 1 ;;
  esac
done
minos=$(otool -l "$qraft" | awk '$1 == "minos" {print $2; exit}')
if [ "$minos" != 12.0 ]; then
  echo "build-macos.sh: the binary needs macOS $minos, not 12.0" >&2
  exit 1
fi
"$qraft" selftest -t 2

dir=$stage/$(archive_name "$target" "$version")
stage_unix "$dir" "$qraft" "$root"
{
  "$qraft" version
  echo
  echo "system: macOS $(sw_vers -productVersion) ($(sw_vers -buildVersion)), $(uname -m)"
  echo "compiler: $(c++ --version | head -n 1)"
  echo "sdk: $(xcrun --show-sdk-version)"
  echo "minimum macOS: $minos"
  echo "libraries:"
  for library in $libraries; do echo "  $library"; done
} > "$dir/BUILDINFO.txt"
pack "$target" "$version" "$stage" "$out" "$root"
ls -l "$out"
