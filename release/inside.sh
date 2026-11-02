#!/usr/bin/env bash
# Runs inside the release image (release/build-release.sh): copies the
# repository to /tmp/src, builds the static qraft, checks it, and writes the
# archive and SHA256SUMS to /out.
set -euo pipefail
umask 022
src=/tmp/src
build=/tmp/build
stage=/tmp/stage
rm -rf "$src" "$build" "$stage"
mkdir -p "$src"
# The committed files, the submodule and the git metadata (for the version).
git -C /repo ls-files --recurse-submodules -z | (cd /repo && xargs -0 cp --parents -t "$src")
cp -a /repo/.git "$src/.git"
version=$(sed -n 's/^project(qraft VERSION \([0-9.]*\).*/\1/p' "$src/CMakeLists.txt")
name=qraft-$version-linux-x86_64

cmake -S "$src" -B "$build" -DCMAKE_BUILD_TYPE=Release -DQRAFT_BUILD_TESTS=OFF \
  "-DCMAKE_CXX_FLAGS=-ffile-prefix-map=$src=. -ffile-prefix-map=$build=."
cmake --build "$build" --target qraft -j "${JOBS:-$(nproc)}"
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
# Modes do not depend on the umask or on the checkout.
find "$stage" -type d -exec chmod 755 {} +
find "$stage" -type f -exec chmod 644 {} +
chmod 755 "$dir/bin/qraft" "$dir/install-check.sh"
{
  "$build/qraft" version
  echo
  echo "image: $(sed -n 's/^FROM //p' "$src/release/Dockerfile")"
  echo "apt: $(sed -n 's/^URIs: //p' /etc/apt/sources.list.d/debian.sources | head -1)"
  echo "packages:"
  dpkg-query -W -f '  ${Package} ${Version}\n' g++ libc6-dev libstdc++-14-dev libjemalloc-dev cmake binutils
} > "$dir/BUILDINFO.txt"

touch -d "@$SOURCE_DATE_EPOCH" $(find "$stage" -print)
tar --sort=name --owner=0 --group=0 --numeric-owner --mtime="@$SOURCE_DATE_EPOCH" \
    -C "$stage" -cf - "$name" | gzip -n -9 > "/out/$name.tar.gz"
(cd /out && sha256sum "$name.tar.gz" > SHA256SUMS)
cat /out/SHA256SUMS
