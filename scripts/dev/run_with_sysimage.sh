#!/usr/bin/env sh

set -eu

BUILD_IF_MISSING=0
STRICT=0

while [ "$#" -gt 0 ]; do
    case "$1" in
        --build-if-missing)
            BUILD_IF_MISSING=1
            shift
            ;;
        --strict)
            STRICT=1
            shift
            ;;
        --help|-h)
            cat <<'EOF'
Usage:
  scripts/dev/run_with_sysimage.sh [--build-if-missing] [--strict] <script-or-julia-args...>

Behavior:
  - if a compatible local sysimage exists, use it
  - otherwise fall back to plain 'julia --project=.'
  - with --build-if-missing, build the local sysimage first when missing
EOF
            exit 0
            ;;
        *)
            break
            ;;
    esac
done

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/../.." && pwd)
SYSIMAGE_PATH="$REPO_ROOT/build/JuliaRelaxTime.so"
META_PATH="$REPO_ROOT/build/JuliaRelaxTime.sysimage.json"

uname_s=$(uname -s 2>/dev/null || printf '')
case "$uname_s" in
    Darwin) SYSIMAGE_PATH="$REPO_ROOT/build/JuliaRelaxTime.dylib" ;;
    Linux) SYSIMAGE_PATH="$REPO_ROOT/build/JuliaRelaxTime.so" ;;
esac

get_julia_version() {
    julia --version | sed -E 's/.* ([0-9]+\.[0-9]+\.[0-9]+).*/\1/'
}

ensure_sysimage() {
    if [ -f "$SYSIMAGE_PATH" ] && [ -f "$META_PATH" ]; then
        return
    fi
    if [ "$BUILD_IF_MISSING" -eq 0 ]; then
        if [ "$STRICT" -eq 1 ]; then
            printf '%s\n' "Sysimage or metadata missing. Re-run with --build-if-missing." >&2
            exit 1
        fi
        return
    fi
    julia --project="$REPO_ROOT" "$REPO_ROOT/scripts/dev/build_sysimage.jl"
}

json_value() {
    key=$1
    sed -n -E "s/.*\"$key\"[[:space:]]*:[[:space:]]*\"([^\"]+)\".*/\1/p" "$META_PATH" | head -n 1
}

ensure_sysimage

USE_SYSIMAGE=0
if [ -f "$SYSIMAGE_PATH" ] && [ -f "$META_PATH" ]; then
    META_VERSION=$(json_value julia_version)
    CURRENT_VERSION=$(get_julia_version)
    if [ "$META_VERSION" = "$CURRENT_VERSION" ]; then
        USE_SYSIMAGE=1
    elif [ "$STRICT" -eq 1 ]; then
        printf '%s\n' "Sysimage Julia version $META_VERSION does not match current Julia $CURRENT_VERSION" >&2
        exit 1
    fi
fi

if [ "$USE_SYSIMAGE" -eq 1 ]; then
    exec julia --sysimage="$SYSIMAGE_PATH" --project="$REPO_ROOT" "$@"
else
    exec julia --project="$REPO_ROOT" "$@"
fi
