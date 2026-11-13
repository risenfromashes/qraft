# Builds the .deb and .rpm of a Linux release from its staged archive folder
# (sourced by inside.sh, after names.sh and stage.sh). Both packages hold the
# same files, from the same staged folder, so their qraft is byte for byte the
# binary of the .tar.gz:
#
#   /usr/bin/qraft
#   /usr/share/man/man1/qraft.1.gz
#   /usr/share/bash-completion/completions/qraft
#   /usr/share/fish/vendor_completions.d/qraft.fish
#   /usr/share/qraft/completions/qraft.zsh   (sourced from ~/.zshrc)
#   /usr/share/doc/qraft/                    (license, notices, README, notes, licenses/)
#
# The .deb also holds a Debian changelog and a lintian override for the
# static binary.
#
# The binary is static, so the packages depend on nothing. Times come from
# SOURCE_DATE_EPOCH.

package_summary="species trees from gene trees by quartet regrafting"
package_description="QRAFT infers a species tree from gene trees: it looks for the species tree
that agrees with the most gene-tree quartets, counting every quartet exactly,
and it scales to a million taxa."
package_maintainer="risenfromashes <ashrafurkhan37@gmail.com>"
package_url="https://github.com/risenfromashes/qraft"

# package_root DIR ROOT lays out the installed files of the staged archive
# folder DIR under ROOT.
package_root() {
  local dir=$1 root=$2
  rm -rf "$root"
  install -D -m 755 "$dir/bin/qraft" "$root/usr/bin/qraft"
  install -d -m 755 "$root/usr/share/man/man1"
  gzip -9 -n -c "$dir/man/man1/qraft.1" > "$root/usr/share/man/man1/qraft.1.gz"
  install -D -m 644 "$dir/completions/qraft.bash" \
    "$root/usr/share/bash-completion/completions/qraft"
  install -D -m 644 "$dir/completions/qraft.fish" \
    "$root/usr/share/fish/vendor_completions.d/qraft.fish"
  install -D -m 644 "$dir/completions/qraft.zsh" \
    "$root/usr/share/qraft/completions/qraft.zsh"
  install -d -m 755 "$root/usr/share/doc/qraft/licenses"
  install -m 644 "$dir/LICENSE" "$dir/README.md" "$dir/THIRD_PARTY_NOTICES.md" \
    "$dir/RELEASE_NOTES.md" "$dir/BUILDINFO.txt" "$root/usr/share/doc/qraft/"
  install -m 644 "$dir"/licenses/* "$root/usr/share/doc/qraft/licenses/"
  {
    # The notice of the README's License section.
    grep -m 1 '^Copyright (C)' "$dir/README.md"
    echo
    echo "QRAFT is free software under the GNU Affero General Public License,"
    echo "version 3 or later: /usr/share/doc/qraft/LICENSE."
    echo
    echo "The binary also contains third-party code; see"
    echo "/usr/share/doc/qraft/THIRD_PARTY_NOTICES.md and the license texts in"
    echo "/usr/share/doc/qraft/licenses/."
  } > "$root/usr/share/doc/qraft/copyright"
  chmod 644 "$root/usr/share/doc/qraft/copyright"
  find "$root" -exec touch -h -d "@$SOURCE_DATE_EPOCH" {} +
}

# make_deb DIR TARGET VERSION OUT writes OUT/<archive name>.deb.
make_deb() {
  local dir=$1 target=$2 version=$3 out=$4 root=/tmp/pkg-deb arch
  case $target in
    linux-x86_64) arch=amd64 ;;
    linux-aarch64) arch=arm64 ;;
  esac
  package_root "$dir" "$root"
  # The binary is static on purpose: it needs no libraries on any distribution.
  install -d -m 755 "$root/usr/share/lintian/overrides"
  echo "qraft: statically-linked-binary [usr/bin/qraft]" \
    > "$root/usr/share/lintian/overrides/qraft"
  chmod 644 "$root/usr/share/lintian/overrides/qraft"
  # A one-entry Debian changelog; the release notes are RELEASE_NOTES.md.
  {
    echo "qraft ($version) unstable; urgency=medium"
    echo
    echo "  * QRAFT $version; see /usr/share/doc/qraft/RELEASE_NOTES.md."
    echo
    echo " -- $package_maintainer  $(LC_ALL=C date -u -R -d "@$SOURCE_DATE_EPOCH")"
  } | gzip -9 -n > "$root/usr/share/doc/qraft/changelog.gz"
  chmod 644 "$root/usr/share/doc/qraft/changelog.gz"
  install -d -m 755 "$root/DEBIAN"
  {
    echo "Package: qraft"
    echo "Version: $version"
    echo "Architecture: $arch"
    echo "Maintainer: $package_maintainer"
    echo "Installed-Size: $(du -sk --apparent-size "$root/usr" | cut -f1)"
    echo "Section: science"
    echo "Priority: optional"
    echo "Homepage: $package_url"
    echo "Description: $package_summary"
    echo "$package_description" | sed 's/^/ /'
  } > "$root/DEBIAN/control"
  chmod 644 "$root/DEBIAN/control"
  find "$root" -exec touch -h -d "@$SOURCE_DATE_EPOCH" {} +
  dpkg-deb --root-owner-group -Zxz -z9 --build "$root" \
    "$out/$(package_file "$target" "$version" deb)"
}

# make_rpm DIR TARGET VERSION OUT writes OUT/<archive name>.rpm.
make_rpm() {
  local dir=$1 target=$2 version=$3 out=$4 root=/tmp/pkg-rpm top=/tmp/rpmbuild arch
  case $target in
    linux-x86_64) arch=x86_64 ;;
    linux-aarch64) arch=aarch64 ;;
  esac
  package_root "$dir" "$root"
  rm -rf "$top"
  mkdir -p "$top/SPECS" "$top/RPMS"
  cat > "$top/SPECS/qraft.spec" <<EOF
%global debug_package %{nil}
%global __os_install_post %{nil}
%global _build_id_links none
%global _buildhost reproducible
%global use_source_date_epoch_as_buildtime 1
%global clamp_mtime_to_source_date_epoch 1
Name: qraft
Version: $version
Release: 1
Summary: QRAFT: $package_summary
License: AGPL-3.0-or-later
URL: $package_url
AutoReqProv: no

%description
$package_description

%install
mkdir -p %{buildroot}
cp -a $root/usr %{buildroot}/

%files
/usr/bin/qraft
/usr/share/man/man1/qraft.1.gz
/usr/share/bash-completion/completions/qraft
/usr/share/fish/vendor_completions.d/qraft.fish
/usr/share/qraft
/usr/share/doc/qraft
EOF
  rpmbuild -bb --target "$arch" --define "_topdir $top" "$top/SPECS/qraft.spec"
  cp "$top/RPMS/$arch/qraft-$version-1.$arch.rpm" \
    "$out/$(package_file "$target" "$version" rpm)"
}
