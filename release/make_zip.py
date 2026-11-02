#!/usr/bin/env python3
"""Writes a reproducible zip of one folder: entries in name order, every time set
to SOURCE_DATE_EPOCH (UTC, at least 1980), modes 755 for folders and executables
and 644 otherwise, deflated at level 9. Usage: make_zip.py ARCHIVE PARENT NAME,
which zips PARENT/NAME as NAME/..."""
import os
import sys
import time
import zipfile

archive, parent, name = sys.argv[1:4]
stamp = time.gmtime(max(int(os.environ.get("SOURCE_DATE_EPOCH", "315532800")), 315532800))[:6]
paths = []
for root, dirs, files in os.walk(os.path.join(parent, name)):
    dirs.sort()
    rel = os.path.relpath(root, parent)
    paths.append((rel + "/", None))
    paths.extend((os.path.join(rel, f), os.path.join(root, f)) for f in sorted(files))
with zipfile.ZipFile(archive, "w") as z:
    for arcname, source in sorted(paths):
        info = zipfile.ZipInfo(arcname.replace(os.sep, "/"), stamp)
        info.create_system = 3
        if source is None:
            info.external_attr = (0o40755 << 16) | 0x10
            z.writestr(info, b"")
            continue
        executable = arcname.endswith((".exe", ".cmd")) or os.access(source, os.X_OK)
        info.external_attr = (0o100755 if executable else 0o100644) << 16
        info.compress_type = zipfile.ZIP_DEFLATED
        with open(source, "rb") as f:
            z.writestr(info, f.read(), compresslevel=9)
