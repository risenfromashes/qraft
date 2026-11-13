# The build environment of the static release binaries: Debian 13 (trixie) at
# a fixed image digest (an index of every architecture's image), with GCC,
# CMake, jemalloc's static library and the .deb and .rpm tools for Linux, and
# on x86-64 also MinGW-w64 GCC for Windows and Wine to run the Windows
# selftest, all from the Debian snapshot of that image's date (debian.sources).
FROM debian:trixie@sha256:2dc419c6d7fe3dad0a2f9fd469c34cc90b0f9ba642d8fde837fe3d7d34a0e7f1
COPY debian.sources /etc/apt/sources.list.d/debian.sources
RUN windows= \
 && if [ "$(dpkg --print-architecture)" = amd64 ]; then \
      windows="g++-mingw-w64-x86-64-posix binutils-mingw-w64-x86-64 wine wine64"; \
    fi \
 && apt-get -o Acquire::Check-Valid-Until=false update \
 && apt-get install -y --no-install-recommends \
      g++ cmake make git python3 libjemalloc-dev binutils file ca-certificates \
      rpm \
      $windows \
 && rm -rf /var/lib/apt/lists/* \
 && git config --system --add safe.directory '*'
