#!/usr/bin/env sh
# Build the documentation (HTML) so you can test that docs work.
# Output: docs/build/. Open docs/build/index.html in a browser to view.
# Usage: ./scripts/build_docs.sh

set -eu
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

if [ ! -f "docs/Project.toml" ]; then
  echo "Error: docs/Project.toml not found." >&2
  exit 1
fi

julia --project=docs -e 'using Pkg; Pkg.instantiate()'
cd docs && julia --project=. -e 'include("make.jl")'
