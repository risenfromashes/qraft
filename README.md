# QRAFT

QRAFT (Quartet Regrafting) infers a species tree from gene trees by
maximizing quartet support, with every quartet counted exactly.

Build with CMake 3.20 or newer and a C++20 compiler:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
build/qraft -i genes.tre -o species.tre
```

Copyright (C) 2026 Ashrafur Rahman, Anik Saha and Md. Shamsuzzoha Bayzid.
QRAFT is free software under the GNU Affero General Public License, version 3
or later; see `LICENSE`.
