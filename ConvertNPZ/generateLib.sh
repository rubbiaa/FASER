#!/usr/bin/env bash
# Convenience wrapper for FASER's existing PyROOT dictionary target.
# Usage: ./generateLib.sh [/path/to/FASER] [-DGENFIT_ROOT=... ...]
set -euo pipefail

faser_source_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
if (( $# > 0 )) && [[ "$1" != -* ]]; then
    faser_source_dir="$1"
    shift
fi

cmake -S "$faser_source_dir" -B "$faser_source_dir/build" "$@"
cmake --build "$faser_source_dir/build" --target CoreUtilsDict \
    --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-8}"
