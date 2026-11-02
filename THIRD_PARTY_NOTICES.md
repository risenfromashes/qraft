# Third-party notices

QRAFT itself is under the GNU Affero General Public License, version 3 or
later (`LICENSE`). The source tree and the release binary also contain the
following third-party code and data.

## Taskflow (source, as the git submodule `external/taskflow`; in the binary)

Used for parallel work. The build overlays one file of it (the worker lookup)
in the build tree; the submodule itself is unchanged.

```
TASKFLOW MIT LICENSE 

Copyright (c) 2018-2026 Dr. Tsung-Wei Huang

The University of Wisconsin at Madison

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## jemalloc (in the release binary, linked statically)

The memory allocator of the release binary. It is not part of the source tree.

```
Copyright: 2002-present Jason Evans <jasone@canonware.com>
          2007-2012 Mozilla Foundation
          2009-present Facebook, Inc.
License: BSD-2-Clause
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice(s),
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice(s),
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER(S) ``AS IS'' AND ANY EXPRESS
OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO
EVENT SHALL THE COPYRIGHT HOLDER(S) BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

## GNU C Library, libstdc++ and libgcc (in the release binary, linked statically)

The release binary is statically linked with the GNU C Library (GNU LGPL,
version 2.1 or later) and with the GCC runtime libraries libstdc++ and libgcc
(GNU GPL, version 3 or later, with the GCC Runtime Library Exception, which
allows their use in binaries under any license). The release archive's
`BUILDINFO.txt` names the exact package versions; their sources are available
from Debian (https://sources.debian.org/) and from https://sourceware.org/glibc/
and https://gcc.gnu.org/. As the LGPL requires, the complete source of QRAFT
is available, so the binary can be relinked with a modified C library by
building QRAFT from source.

The release archive's `licenses/` folder holds the full texts:
`GPL-3.0.txt`, `LGPL-2.1.txt`, `gcc-runtime-libraries.txt` (the GCC
copyright file, with the GCC Runtime Library Exception), `glibc.txt`,
`jemalloc.txt` and `taskflow.txt`.

## Test data

`tests/trees/astral2_200taxa_*.tre` are derived from the simulated 200-taxon
gene trees of the ASTRAL-II study (Mirarab and Warnow, Bioinformatics, 2015).

## Example data

`examples/` holds gene trees we simulated with SimPhy (Mallo, de Oliveira
Martins and Posada, Systematic Biology, 2016): 100 taxa, the first 200 of
1,000 genes, branch lengths removed.
