#!/usr/bin/env python3
"""Writes a reproducible release archive of one folder.

Usage: make_archive.py ARCHIVE PARENT NAME

Packs PARENT/NAME as NAME/... into ARCHIVE, a .tar.gz or a .zip by its
extension. Entries come in name order, every time is SOURCE_DATE_EPOCH (UTC,
at least 1980), owners are 0, and modes are 755 for folders, for files in a
bin/ folder and for .sh, .cmd and .exe files, and 644 for the rest, whatever
the checkout's modes were. The tar is gzipped at level 9 with no name or time
in the gzip header; the zip is deflated at level 9. The same files give the
same bytes with the same Python and zlib.
"""
import gzip
import io
import os
import sys
import tarfile
import time
import zipfile


def entries(parent, name):
    """(archive name, source path or None for a folder), sorted by name."""
    found = []
    for root, dirs, files in os.walk(os.path.join(parent, name)):
        dirs.sort()
        rel = os.path.relpath(root, parent).replace(os.sep, "/")
        found.append((rel + "/", None))
        found.extend((rel + "/" + f, os.path.join(root, f)) for f in sorted(files))
    return sorted(found)


def executable(arcname):
    return "/bin/" in arcname or arcname.endswith((".sh", ".cmd", ".exe"))


def write_tar_gz(archive, items, epoch):
    data = io.BytesIO()
    with tarfile.open(fileobj=data, mode="w", format=tarfile.USTAR_FORMAT) as tar:
        for arcname, source in items:
            info = tarfile.TarInfo(arcname.rstrip("/"))
            info.mtime = epoch
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            if source is None:
                info.type = tarfile.DIRTYPE
                info.mode = 0o755
                tar.addfile(info)
                continue
            with open(source, "rb") as f:
                content = f.read()
            info.size = len(content)
            info.mode = 0o755 if executable(arcname) else 0o644
            tar.addfile(info, io.BytesIO(content))
    with open(archive, "wb") as out:
        with gzip.GzipFile(filename="", mode="wb", fileobj=out,
                           compresslevel=9, mtime=0) as gz:
            gz.write(data.getvalue())


def write_zip(archive, items, epoch):
    stamp = time.gmtime(max(epoch, 315532800))[:6]
    with zipfile.ZipFile(archive, "w") as z:
        for arcname, source in items:
            info = zipfile.ZipInfo(arcname, stamp)
            info.create_system = 3
            if source is None:
                info.external_attr = (0o40755 << 16) | 0x10
                z.writestr(info, b"")
                continue
            mode = 0o100755 if executable(arcname) else 0o100644
            info.external_attr = mode << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            with open(source, "rb") as f:
                z.writestr(info, f.read(), compresslevel=9)


def main():
    archive, parent, name = sys.argv[1:4]
    epoch = int(os.environ.get("SOURCE_DATE_EPOCH", "315532800"))
    items = entries(parent, name)
    if archive.endswith(".tar.gz"):
        write_tar_gz(archive, items, epoch)
    elif archive.endswith(".zip"):
        write_zip(archive, items, epoch)
    else:
        sys.exit("make_archive.py: ARCHIVE must end in .tar.gz or .zip")


if __name__ == "__main__":
    main()
