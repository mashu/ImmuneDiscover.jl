#!/usr/bin/env sh
set -eu

# Run immunediscover without PackageCompiler
# Usage: ./run-immunediscover.sh <args...>
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR"

# Prefer bundled Julia if available, otherwise use system julia
# Prefer system Julia first to avoid loading a baked sysimage; fallback to bundled
if command -v julia >/dev/null 2>&1; then
  JULIA_BIN="$(command -v julia)"
elif [ -x "$PROJECT_DIR/immunediscover/bin/julia" ]; then
  JULIA_BIN="$PROJECT_DIR/immunediscover/bin/julia"
else
  echo "Error: julia not found" >&2
  exit 1
fi

# Force-load source module to avoid any precompiled sysimage overriding local changes
EXPR="include(\"$PROJECT_DIR/src/immunediscover.jl\"); exit(immunediscover.julia_main())"
exec "$JULIA_BIN" --project="$PROJECT_DIR" -e "$EXPR" -- "$@"


