#!/bin/sh
# Installs this unpacked QRAFT release: the qraft binary, its man page and
# shell completions, under a prefix (default ~/.local; /usr/local needs root).
#
#   ./install.sh [--prefix DIR]
#   ./install.sh [--prefix DIR] --uninstall
set -eu
here=$(cd "$(dirname "$0")" && pwd)
prefix=$HOME/.local
uninstall=0
while [ $# -gt 0 ]; do
  case $1 in
    --prefix) [ $# -ge 2 ] || { echo "install.sh: --prefix needs a folder" >&2; exit 2; }
              prefix=$2; shift 2 ;;
    --prefix=*) prefix=${1#--prefix=}; shift ;;
    --uninstall) uninstall=1; shift ;;
    -h | --help) sed -n '2,6p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "install.sh: unknown argument '$1' (see --help)" >&2; exit 2 ;;
  esac
done

# Each line: the file in this folder, then where it goes under the prefix.
files="bin/qraft bin/qraft
man/man1/qraft.1 share/man/man1/qraft.1
completions/qraft.bash share/bash-completion/completions/qraft
completions/qraft.fish share/fish/vendor_completions.d/qraft.fish
completions/qraft.zsh share/qraft/completions/qraft.zsh"

echo "$files" | while read -r from to; do
  dest=$prefix/$to
  if [ "$uninstall" = 1 ]; then
    rm -f "$dest"
  else
    mkdir -p "$(dirname "$dest")"
    mode=644
    [ "$from" = bin/qraft ] && mode=755
    # Copy beside the old file, then rename over it: a running qraft keeps
    # its old file, where copying onto it would fail ("text file busy").
    tmp=$dest.install-$$
    if ! { cp "$here/$from" "$tmp" && chmod "$mode" "$tmp" && mv -f "$tmp" "$dest"; }; then
      rm -f "$tmp"
      echo "install.sh: cannot install $dest" >&2
      exit 1
    fi
  fi
done

if [ "$uninstall" = 1 ]; then
  rmdir "$prefix/share/qraft/completions" "$prefix/share/qraft" 2>/dev/null || true
  echo "removed qraft from $prefix"
  exit 0
fi
echo "installed $("$prefix/bin/qraft" --version) in $prefix"
case :$PATH: in
  *:"$prefix/bin":*) ;;
  *) echo "add $prefix/bin to your PATH, for example in ~/.bashrc:"
     echo "  export PATH=\"$prefix/bin:\$PATH\"" ;;
esac
echo "zsh completion: add 'source $prefix/share/qraft/completions/qraft.zsh' to ~/.zshrc"
