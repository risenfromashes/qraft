#!/usr/bin/env bash
# Runs inside the release image (release/build-release.sh): copies the
# repository to /tmp/src, builds the static Linux and Windows binaries, checks
# them, and writes the two archives and SHA256SUMS to /out.
set -euo pipefail
umask 022
src=/tmp/src
build=/tmp/build
wbuild=/tmp/build-windows
stage=/tmp/stage
rm -rf "$src" "$build" "$wbuild" "$stage"
mkdir -p "$src"
# The committed files, the submodule and the git metadata (for the version).
git -C /repo ls-files --recurse-submodules -z | (cd /repo && xargs -0 cp --parents -t "$src")
cp -a /repo/.git "$src/.git"
version=$(sed -n 's/^project(qraft VERSION \([0-9.]*\).*/\1/p' "$src/CMakeLists.txt")
jobs=${JOBS:-$(nproc)}
image=$(sed -n 's/^FROM //p' "$src/release/Dockerfile")
apt_source=$(sed -n 's/^URIs: //p' /etc/apt/sources.list.d/debian.sources | head -1)

# Linux x86-64.
name=qraft-$version-linux-x86_64
cmake -S "$src" -B "$build" -DCMAKE_BUILD_TYPE=Release -DQRAFT_BUILD_TESTS=OFF \
  "-DCMAKE_CXX_FLAGS=-ffile-prefix-map=$src=. -ffile-prefix-map=$build=."
cmake --build "$build" --target qraft -j "$jobs"
strip --strip-all "$build/qraft"
"$build/qraft" selftest -t 2

dir=$stage/$name
mkdir -p "$dir/bin" "$dir/man/man1" "$dir/completions" "$dir/examples"
cp "$build/qraft" "$dir/bin/qraft"
"$build/qraft" help man > "$dir/man/man1/qraft.1"
"$build/qraft" completion bash > "$dir/completions/qraft.bash"
"$build/qraft" completion zsh > "$dir/completions/qraft.zsh"
"$build/qraft" completion fish > "$dir/completions/qraft.fish"
cp "$src/LICENSE" "$src/README.md" "$src/THIRD_PARTY_NOTICES.md" "$src/CITATION.cff" \
   "$src/release/RELEASE_NOTES.md" "$dir/"
cp "$src/release/install-check.sh" "$dir/"
cp "$src"/examples/* "$dir/examples/"
# The license texts of the code linked into the binary.
mkdir -p "$dir/licenses"
cp /usr/share/common-licenses/GPL-3 "$dir/licenses/GPL-3.0.txt"
cp /usr/share/common-licenses/LGPL-2.1 "$dir/licenses/LGPL-2.1.txt"
cp /usr/share/doc/gcc-14-base/copyright "$dir/licenses/gcc-runtime-libraries.txt"
cp /usr/share/doc/libc6/copyright "$dir/licenses/glibc.txt"
cp /usr/share/doc/libjemalloc2/copyright "$dir/licenses/jemalloc.txt"
cp "$src/external/taskflow/LICENSE" "$dir/licenses/taskflow.txt"
{
  "$build/qraft" version
  echo
  echo "image: $image"
  echo "apt: $apt_source"
  echo "packages:"
  dpkg-query -W -f '  ${Package} ${Version}\n' g++ libc6-dev libstdc++-14-dev libjemalloc-dev cmake binutils
} > "$dir/BUILDINFO.txt"

# Windows x86-64 (MinGW-w64, static, system allocator).
wname=qraft-$version-windows-x86_64
cmake -S "$src" -B "$wbuild" -DCMAKE_SYSTEM_NAME=Windows -DCMAKE_SYSTEM_PROCESSOR=x86_64 \
  -DCMAKE_CXX_COMPILER=x86_64-w64-mingw32-g++-posix \
  -DCMAKE_BUILD_TYPE=Release -DQRAFT_BUILD_TESTS=OFF -DQRAFT_USE_JEMALLOC=OFF \
  "-DCMAKE_CXX_FLAGS=-ffile-prefix-map=$src=. -ffile-prefix-map=$wbuild=."
cmake --build "$wbuild" --target qraft -j "$jobs"
x86_64-w64-mingw32-strip --strip-all "$wbuild/qraft.exe"
mkdir -p /tmp/wine-home
wine() { HOME=/tmp/wine-home WINEPREFIX=/tmp/wine-home/prefix WINEDEBUG=-all command wine "$@" 2>/dev/null; }
wine "$wbuild/qraft.exe" selftest -t 2

wdir=$stage/$wname
mkdir -p "$wdir/examples" "$wdir/licenses"
cp "$wbuild/qraft.exe" "$wdir/qraft.exe"
cp "$src/LICENSE" "$src/README.md" "$src/THIRD_PARTY_NOTICES.md" "$src/CITATION.cff" \
   "$src/release/RELEASE_NOTES.md" "$src/release/install-check.cmd" "$wdir/"
cp "$src"/examples/* "$wdir/examples/"
cp /usr/share/common-licenses/GPL-3 "$wdir/licenses/GPL-3.0.txt"
cp /usr/share/doc/gcc-mingw-w64-base/copyright "$wdir/licenses/gcc-runtime-libraries.txt"
cp /usr/share/doc/mingw-w64-x86-64-dev/copyright "$wdir/licenses/mingw-w64-runtime.txt"
cp "$src/external/taskflow/LICENSE" "$wdir/licenses/taskflow.txt"
{
  wine "$wbuild/qraft.exe" version | tr -d '\r'
  echo
  echo "image: $image"
  echo "apt: $apt_source"
  echo "packages:"
  dpkg-query -W -f '  ${Package} ${Version}\n' g++-mingw-w64-x86-64-posix mingw-w64-x86-64-dev \
    binutils-mingw-w64-x86-64 cmake
} > "$wdir/BUILDINFO.txt"

# Modes do not depend on the umask or on the checkout.
find "$stage" -type d -exec chmod 755 {} +
find "$stage" -type f -exec chmod 644 {} +
chmod 755 "$dir/bin/qraft" "$dir/install-check.sh" "$wdir/qraft.exe"
touch -d "@$SOURCE_DATE_EPOCH" $(find "$stage" -print)
tar --sort=name --owner=0 --group=0 --numeric-owner --mtime="@$SOURCE_DATE_EPOCH" \
    -C "$stage" -cf - "$name" | gzip -n -9 > "/out/$name.tar.gz"
python3 "$src/release/make_zip.py" "/out/$wname.zip" "$stage" "$wname"
(cd /out && sha256sum "$name.tar.gz" "$wname.zip" > SHA256SUMS)
cat /out/SHA256SUMS
