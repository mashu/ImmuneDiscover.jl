#!/usr/bin/env sh
set -eu
# Run immunediscover from source (no package build). Usage: ./scripts/run.sh [args...]
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
if command -v julia >/dev/null 2>&1; then
  JULIA_BIN="$(command -v julia)"
else
  echo "Error: julia not found" >&2
  exit 1
fi
# `using` loads the precompiled package image (@compile_workload in src/immunediscover.jl).
# `include()` bypasses that cache and re-JITs the whole module on every invocation.
exec "$JULIA_BIN" --project="$PROJECT_DIR" -e 'using immunediscover; exit(immunediscover.julia_main())' -- "$@"
