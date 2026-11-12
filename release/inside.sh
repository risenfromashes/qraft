#!/usr/bin/env bash
# Runs inside the release image (release/build-release.sh): copies the
# repository to /tmp/src, then builds, checks and packs the static binary of
# each TARGET given, and writes the archives to /out.
#
#   inside.sh [TARGET...]
#
# TARGET is linux-x86_64 or windows-x86_64 on an x86-64 host and
# linux-aarch64 on an arm64 host; without one, the host's targets are built.
set -euo pipefail
umask 022
src=/tmp/src
stage=/tmp/stage
rm -rf "$src" "$stage" /tmp/build-*
mkdir -p "$src" "$stage"
# The committed files, the submodule and the git metadata (for the version).
git -C /repo ls-files --recurse-submodules -z | (cd /repo && xargs -0 cp --parents -t "$src")
cp -a /repo/.git "$src/.git"
. "$src/release/names.sh"
. "$src/release/stage.sh"
version=$(release_version "$src")
jobs=${JOBS:-$(nproc)}
image=$(sed -n 's/^FROM //p' "$src/release/Dockerfile")
apt_source=$(sed -n 's/^URIs: //p' /etc/apt/sources.list.d/debian.sources | head -1)

# build_info PACKAGE... prints the image, the apt snapshot and the versions of
# the named packages, for BUILDINFO.txt.
build_info() {
  echo
  echo "image: $image"
  echo "apt: $apt_source"
  echo "packages:"
  dpkg-query -W -f '  ${Package} ${Version}\n' "$@"
}

# Linux, natively on the host's architecture: static, with jemalloc.
build_linux() {
  local target=$1 build=/tmp/build-$1 dir
  case $target/$(uname -m) in
    linux-x86_64/x86_64 | linux-aarch64/aarch64) ;;
    *) echo "inside.sh: $target cannot be built on $(uname -m)" >&2; exit 2 ;;
  esac
  cmake -S "$src" -B "$build" -DCMAKE_BUILD_TYPE=Release -DQRAFT_BUILD_TESTS=OFF \
    "-DCMAKE_CXX_FLAGS=-ffile-prefix-map=$src=. -ffile-prefix-map=$build=."
  cmake --build "$build" --target qraft -j "$jobs"
  strip --strip-all "$build/qraft"
  "$build/qraft" selftest -t 2

  dir=$stage/$(archive_name "$target" "$version")
  stage_unix "$dir" "$build/qraft" "$src"
  # The license texts of the code linked into the binary.
  cp /usr/share/common-licenses/GPL-3 "$dir/licenses/GPL-3.0.txt"
  cp /usr/share/common-licenses/LGPL-2.1 "$dir/licenses/LGPL-2.1.txt"
  cp /usr/share/doc/gcc-14-base/copyright "$dir/licenses/gcc-runtime-libraries.txt"
  cp /usr/share/doc/libc6/copyright "$dir/licenses/glibc.txt"
  cp /usr/share/doc/libjemalloc2/copyright "$dir/licenses/jemalloc.txt"
  {
    "$build/qraft" version
    build_info g++ libc6-dev libstdc++-14-dev libjemalloc-dev cmake binutils
  } > "$dir/BUILDINFO.txt"
  pack "$target" "$version" "$stage" /out "$src"
}

# Windows x86-64: cross-compiled with MinGW-w64, static, with the system
# allocator, and checked under Wine.
build_windows() {
  local target=windows-x86_64 build=/tmp/build-windows dir
  cmake -S "$src" -B "$build" -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_SYSTEM_PROCESSOR=x86_64 \
    -DCMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++-posix \
    -DCMAKE_BUILD_TYPE=Release -DQRAFT_BUILD_TESTS=OFF -DQRAFT_USE_JEMALLOC=OFF \
    "-DCMAKE_CXX_FLAGS=-ffile-prefix-map=$src=. -ffile-prefix-map=$build=."
  cmake --build "$build" --target qraft -j "$jobs"
  x86_64-w64-mingw32-strip --strip-all "$build/qraft.exe"
  mkdir -p /tmp/wine-home
  wine() { HOME=/tmp/wine-home WINEPREFIX=/tmp/wine-home/prefix WINEDEBUG=-all command wine "$@" 2>/dev/null; }
  wine "$build/qraft.exe" selftest -t 2

  dir=$stage/$(archive_name "$target" "$version")
  mkdir -p "$dir"
  cp "$build/qraft.exe" "$dir/qraft.exe"
  stage_common "$dir" "$src"
  cp "$src/release/install-check.cmd" "$dir/"
  cp /usr/share/common-licenses/GPL-3 "$dir/licenses/GPL-3.0.txt"
  cp /usr/share/doc/gcc-mingw-w64-base/copyright "$dir/licenses/gcc-runtime-libraries.txt"
  cp /usr/share/doc/mingw-w64-x86-64-dev/copyright "$dir/licenses/mingw-w64-runtime.txt"
  {
    wine "$build/qraft.exe" version | tr -d '\r'
    build_info g++-mingw-w64-x86-64-posix mingw-w64-x86-64-dev binutils-mingw-w64-x86-64 cmake
  } > "$dir/BUILDINFO.txt"
  pack "$target" "$version" "$stage" /out "$src"
}

if [ $# -eq 0 ]; then
  case $(uname -m) in
    x86_64) set -- linux-x86_64 windows-x86_64 ;;
    aarch64) set -- linux-aarch64 ;;
    *) echo "inside.sh: no release targets for $(uname -m)" >&2; exit 2 ;;
  esac
fi
for target in "$@"; do
  case $target in
    linux-x86_64 | linux-aarch64) build_linux "$target" ;;
    windows-x86_64) build_windows ;;
    *) echo "inside.sh: unknown target $target" >&2; exit 2 ;;
  esac
done
