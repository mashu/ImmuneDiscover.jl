#!/usr/bin/env sh
# Build immunediscover standalone binary for the current platform.
# Requires: Julia (1.9+), build/ project with PackageCompiler.
# Usage: ./scripts/build_binary.sh [output_dir]
#   If output_dir is omitted, uses build/immunediscover_app/

set -eu
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

if [ ! -f "build/Project.toml" ]; then
  echo "Error: build/Project.toml not found. Ensure build environment exists." >&2
  exit 1
fi

julia --project=build -e 'using Pkg; Pkg.instantiate()'
julia --project=build scripts/build_binary.jl "$@"
