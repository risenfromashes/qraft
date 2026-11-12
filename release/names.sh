# The name of every release archive, set in one place. The release scripts
# source this file; the workflow and SHA256SUMS take whatever files the scripts
# write, so a change here renames everything (the README's download links
# name the files too).
#
# {target} is linux-x86_64, linux-aarch64, macos-arm64 or windows-x86_64 and
# {version} the version in CMakeLists.txt. The archive unpacks into a folder
# of the same name.
qraft_archive_pattern='qraft-{target}'

# archive_name TARGET VERSION prints the archive's name without its extension,
# which is also the name of the folder inside it.
archive_name() {
  local name=${qraft_archive_pattern//\{target\}/$1}
  printf '%s\n' "${name//\{version\}/$2}"
}

# archive_file TARGET VERSION prints the archive's file name: .zip for
# Windows, .tar.gz otherwise.
archive_file() {
  case $1 in
    windows-*) printf '%s.zip\n' "$(archive_name "$1" "$2")" ;;
    *) printf '%s.tar.gz\n' "$(archive_name "$1" "$2")" ;;
  esac
}
