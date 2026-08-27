# Map argv to cached ArgParse help files. Names match Cli.argparse_named_nodes:
# root, <group>, <group>-<sub>. Caller sets HELP_DIR.

is_help_flag() {
    [ "$1" = "--help" ] || [ "$1" = "-h" ]
}

is_help_request() {
    [ $# -eq 0 ] && return 0
    is_help_flag "$1" && return 0
    [ $# -ge 2 ] && is_help_flag "$2" && return 0
    [ $# -ge 3 ] && is_help_flag "$3" && return 0
    return 1
}

help_page_path() {
    if [ $# -eq 0 ] || { [ $# -eq 1 ] && is_help_flag "$1"; }; then
        printf '%s\n' "$HELP_DIR/root.txt"
        return
    fi
    if [ $# -eq 2 ] && is_help_flag "$2"; then
        printf '%s\n' "$HELP_DIR/$1.txt"
        return
    fi
    if [ $# -eq 3 ] && is_help_flag "$3"; then
        printf '%s\n' "$HELP_DIR/$1-$2.txt"
        return
    fi
}

print_cached_help() {
    is_help_request "$@" || return 1
    page=$(help_page_path "$@")
    [ -n "$page" ] && [ -f "$page" ] || return 1
    cat "$page"
    return 0
}
