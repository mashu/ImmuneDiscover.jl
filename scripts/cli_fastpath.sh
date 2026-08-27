# Instant --help / --version without starting Julia.
# Source after setting HELP_DIR. Caller must define print_cli_version.
# Usage: try_fast_help_or_version "$@" && exit 0

wants_static_help() {
    n=$#
    if [ "$n" -eq 0 ]; then
        return 0
    fi
    case "$1" in --help|-h) return 0 ;; esac
    if [ "$n" -ge 2 ]; then
        case "$2" in --help|-h) return 0 ;; esac
    fi
    if [ "$n" -ge 3 ]; then
        case "$3" in --help|-h) return 0 ;; esac
    fi
    return 1
}

try_fast_help_or_version() {
    helpfile=""
    n=$#

    if [ "$n" -eq 0 ]; then
        helpfile="$HELP_DIR/root.txt"
    elif [ "$n" -eq 1 ]; then
        case "$1" in
            --version|-V) print_cli_version; return 0 ;;
            --help|-h) helpfile="$HELP_DIR/root.txt" ;;
        esac
    elif [ "$n" -eq 2 ]; then
        case "$1" in
            --version|-V) print_cli_version; return 0 ;;
            --help|-h) helpfile="$HELP_DIR/root.txt" ;;
        esac
        case "$2" in
            --version|-V) print_cli_version; return 0 ;;
            --help|-h) helpfile="$HELP_DIR/$1.txt" ;;
        esac
    elif [ "$n" -eq 3 ]; then
        case "$1" in
            --version|-V) print_cli_version; return 0 ;;
            --help|-h) helpfile="$HELP_DIR/root.txt" ;;
        esac
        case "$3" in
            --help|-h) helpfile="$HELP_DIR/$1-$2.txt" ;;
        esac
    fi

    if [ -n "${helpfile:-}" ] && [ -f "$helpfile" ]; then
        cat "$helpfile"
        return 0
    fi
    return 1
}
