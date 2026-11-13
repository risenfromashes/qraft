#!/usr/bin/env bash
# Installs the .deb and .rpm of one Linux target in clean containers (Debian
# 13, Ubuntu 24.04, Fedora, Rocky Linux 9) of the host's architecture and
# checks each install: /usr/bin/qraft is byte for byte the .tar.gz's binary,
# its selftest passes, it infers the example's expected tree, the package
# lists the man page, completions and license (minimal images such as
# Ubuntu's and Fedora's leave documentation out on disk), and removing the
# package removes it.
#
#   release/test-packages.sh DIR TARGET
#
# DIR holds the release files (as written by build-release.sh); needs Docker.
set -euo pipefail
dir=$(cd "$1" && pwd)
target=$2
root=$(git -C "$(dirname "$0")/.." rev-parse --show-toplevel)
. "$root/release/names.sh"
. "$root/release/stage.sh"
version=$(release_version "$root")
name=$(archive_name "$target" "$version")
want=$(tar -xzOf "$dir/$(archive_file "$target" "$version")" "$name/bin/qraft" |
  sha256sum | cut -d ' ' -f 1)
deb=$(package_file "$target" "$version" deb)
rpm=$(package_file "$target" "$version" rpm)

# The checks, run inside each container after the install.
check=$(mktemp)
trap 'rm -f "$check"' EXIT
cat > "$check" <<'EOF'
set -eu
test "$(sha256sum /usr/bin/qraft | cut -d ' ' -f 1)" = "$WANT"
qraft --version
qraft selftest -t 2 > /dev/null
qraft -i /examples/genes-100taxa.tre -o /tmp/out.tre -q
test "$(sha256sum < /tmp/out.tre)" = "$(sha256sum < /examples/expected-100taxa.tre)"
for f in /usr/bin/qraft /usr/share/man/man1/qraft.1.gz \
         /usr/share/bash-completion/completions/qraft \
         /usr/share/fish/vendor_completions.d/qraft.fish \
         /usr/share/qraft/completions/qraft.zsh /usr/share/doc/qraft/LICENSE \
         /usr/share/doc/qraft/licenses/taskflow.txt; do
  grep -qx "$f" /tmp/package-files
done
test -s /usr/share/bash-completion/completions/qraft
EOF
chmod 644 "$check"

# run IMAGE INSTALL LIST REMOVE installs, lists the package's files, checks
# and removes in a fresh container.
run() {
  docker run --rm -e WANT="$want" -e DEBIAN_FRONTEND=noninteractive \
    -v "$dir:/dist:ro" -v "$root/examples:/examples:ro" -v "$check:/check.sh:ro" \
    "$1" sh -c "$2 && $3 > /tmp/package-files && sh /check.sh && $4 && ! test -e /usr/bin/qraft"
  echo "$1: installed, checked and removed"
}
for image in debian:trixie ubuntu:24.04; do
  run "$image" "apt-get install -y -qq /dist/$deb > /dev/null" "dpkg -L qraft" \
    "apt-get remove -y -qq qraft > /dev/null"
done
for image in fedora:latest rockylinux:9; do
  run "$image" "dnf install -y -q --disablerepo='*' /dist/$rpm" "rpm -ql qraft" \
    "dnf remove -y -q qraft"
done
