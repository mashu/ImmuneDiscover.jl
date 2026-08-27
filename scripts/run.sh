#!/usr/bin/env sh
set -eu
# Run immunediscover from source. Usage: ./scripts/run.sh [args...]
SCRIPT_DIR="$(CDPATH= cd -- "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)"
HELP_DIR="$PROJECT_DIR/build/help"
SYSIMAGE="$PROJECT_DIR/build/immunediscover.so"

# shellcheck source=cli_fastpath.sh
. "$SCRIPT_DIR/cli_fastpath.sh"

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

help_cache_stale() {
    [ -f "$HELP_DIR/root.txt" ] || return 0
    for f in "$PROJECT_DIR"/src/cmd/*.jl "$PROJECT_DIR"/Project.toml; do
        [ -f "$f" ] || continue
        [ "$f" -nt "$HELP_DIR/root.txt" ] && return 0
    done
    return 1
}

find_julia() {
    if command -v julia >/dev/null 2>&1; then
        JULIA_BIN="$(command -v julia)"
        return
    fi
    echo "Error: julia not found" >&2
    exit 1
}

invoke_julia() {
    if [ -f "$SYSIMAGE" ]; then
        "$JULIA_BIN" --sysimage="$SYSIMAGE" --startup-file=no --quiet --project="$PROJECT_DIR" "$@"
    else
        "$JULIA_BIN" --startup-file=no --quiet --project="$PROJECT_DIR" "$@"
    fi
}

case "${1:-}" in
    --version|-V) print_cli_version; exit 0 ;;
esac

if is_help_request "$@"; then
    if help_cache_stale; then
        find_julia
        invoke_julia -e 'using immunediscover; immunediscover.Cli.ensure_help_pages!()' >/dev/null
    fi
    print_cached_help "$@" && exit 0
fi

find_julia
if [ -f "$SYSIMAGE" ]; then
    exec "$JULIA_BIN" --sysimage="$SYSIMAGE" --startup-file=no --quiet --project="$PROJECT_DIR" \
        -e 'using immunediscover; exit(immunediscover.julia_main())' -- "$@"
fi
exec "$JULIA_BIN" --startup-file=no --quiet --project="$PROJECT_DIR" \
    -e 'using immunediscover; exit(immunediscover.julia_main())' -- "$@"
