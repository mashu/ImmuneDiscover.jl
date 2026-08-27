#!/usr/bin/env sh
set -eu
# Run immunediscover from source (no package build). Usage: ./scripts/run.sh [args...]
SCRIPT_DIR="$(CDPATH= cd -- "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)"
HELP_DIR="$PROJECT_DIR/build/help"
SYSIMAGE="$PROJECT_DIR/build/immunediscover.so"

print_cli_version() {
    ver=""
    if [ -f "$PROJECT_DIR/Project.toml" ]; then
        ver=$(sed -n 's/^version = "\(.*\)"/\1/p' "$PROJECT_DIR/Project.toml" | head -n 1)
    fi
    hash="unknown"
    if command -v git >/dev/null 2>&1; then
        hash=$(git -C "$PROJECT_DIR" rev-parse HEAD 2>/dev/null || echo unknown)
    fi
    printf '%s (git %s)\n' "$ver" "$hash"
}

# shellcheck source=cli_fastpath.sh
. "$SCRIPT_DIR/cli_fastpath.sh"

# --version never needs ArgParse pages
case "${1:-}" in
    --version|-V) print_cli_version; exit 0 ;;
esac

help_pages_stale() {
    [ -f "$HELP_DIR/root.txt" ] || return 0
    for f in "$PROJECT_DIR"/src/cmd/*.jl "$PROJECT_DIR"/Project.toml; do
        [ -f "$f" ] || continue
        if [ "$f" -nt "$HELP_DIR/root.txt" ]; then
            return 0
        fi
    done
    return 1
}

if ! help_pages_stale && try_fast_help_or_version "$@"; then
    exit 0
fi

if command -v julia >/dev/null 2>&1; then
    JULIA_BIN="$(command -v julia)"
else
    echo "Error: julia not found" >&2
    exit 1
fi

run_julia() {
    if [ -f "$SYSIMAGE" ]; then
        "$JULIA_BIN" --sysimage="$SYSIMAGE" --startup-file=no --quiet --project="$PROJECT_DIR" "$@"
    else
        "$JULIA_BIN" --startup-file=no --quiet --project="$PROJECT_DIR" "$@"
    fi
}

if wants_static_help "$@" && help_pages_stale; then
    mkdir -p "$HELP_DIR"
    run_julia -e 'using immunediscover; immunediscover.Cli.write_cli_help_pages!(ARGS[1])' -- "$HELP_DIR" >/dev/null
    if try_fast_help_or_version "$@"; then
        exit 0
    fi
fi

# Prefer a local sysimage so CSV/DataFrames/FASTX are already native. Rebuild with
# ./scripts/build_sysimage.sh after changing source or the Manifest.
# `using` loads the precompiled package image unless a sysimage already contains it.
if [ -f "$SYSIMAGE" ]; then
    exec "$JULIA_BIN" --sysimage="$SYSIMAGE" --startup-file=no --quiet --project="$PROJECT_DIR" \
        -e 'using immunediscover; exit(immunediscover.julia_main())' -- "$@"
fi
exec "$JULIA_BIN" --startup-file=no --quiet --project="$PROJECT_DIR" \
    -e 'using immunediscover; exit(immunediscover.julia_main())' -- "$@"
