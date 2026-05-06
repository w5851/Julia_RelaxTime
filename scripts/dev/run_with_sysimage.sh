#!/usr/bin/env sh

set -eu

BUILD_IF_MISSING=0
MISMATCH_POLICY=rebuild

while [ "$#" -gt 0 ]; do
    case "$1" in
        --build-if-missing)
            BUILD_IF_MISSING=1
            shift
            ;;
        --mismatch-policy=*)
            MISMATCH_POLICY=${1#*=}
            shift
            ;;
        --help|-h)
            cat <<'EOF'
Usage:
  scripts/dev/run_with_sysimage.sh [--build-if-missing] [--mismatch-policy=fallback|strict|rebuild] <script-or-julia-args...>

Behavior:
  - fallback: use sysimage only when compatible; otherwise fall back to plain 'julia --project=.'
  - strict: require a compatible sysimage, otherwise exit with error
  - rebuild: build a fresh local sysimage when missing or incompatible
EOF
            exit 0
            ;;
        *)
            break
            ;;
    esac
done

if [ "$BUILD_IF_MISSING" -eq 1 ]; then
    if [ "$MISMATCH_POLICY" = "strict" ]; then
        printf '%s\n' "--build-if-missing cannot be combined with --mismatch-policy=strict" >&2
        exit 1
    fi
    MISMATCH_POLICY=rebuild
fi

case "$MISMATCH_POLICY" in
    fallback|strict|rebuild) ;;
    *)
        printf '%s\n' "Unsupported mismatch policy: $MISMATCH_POLICY" >&2
        exit 1
        ;;
esac

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/../.." && pwd)
SYSIMAGE_PATH="$REPO_ROOT/build/JuliaRelaxTime.so"
META_PATH="$REPO_ROOT/build/JuliaRelaxTime.sysimage.json"

uname_s=$(uname -s 2>/dev/null || printf '')
case "$uname_s" in
    Darwin) SYSIMAGE_PATH="$REPO_ROOT/build/JuliaRelaxTime.dylib" ; PLATFORM_FAMILY=macos ;;
    Linux) SYSIMAGE_PATH="$REPO_ROOT/build/JuliaRelaxTime.so" ; PLATFORM_FAMILY=linux ;;
    MINGW*|MSYS*|CYGWIN*) SYSIMAGE_PATH="$REPO_ROOT/build/JuliaRelaxTime.dll" ; PLATFORM_FAMILY=windows ;;
    *) PLATFORM_FAMILY=unknown ;;
esac

get_julia_version() {
    julia --version | sed -E 's/.* ([0-9]+\.[0-9]+\.[0-9]+).*/\1/'
}

get_platform_arch() {
    uname -m 2>/dev/null | sed -E 's/amd64/x86_64/; s/x64/x86_64/; s/arm64/aarch64/'
}

build_local_sysimage() {
    julia --project="$REPO_ROOT" "$REPO_ROOT/scripts/dev/build_sysimage.jl"
}

get_head_commit() {
    if command -v git >/dev/null 2>&1; then
        git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || true
    fi
}

json_value() {
    key=$1
    sed -n -E "s/.*\"$key\"[[:space:]]*:[[:space:]]*\"([^\"]+)\".*/\1/p" "$META_PATH" | head -n 1
}

handle_missing_sysimage() {
    case "$MISMATCH_POLICY" in
        strict)
            printf '%s\n' "Sysimage or metadata missing." >&2
            exit 1
            ;;
        rebuild)
            build_local_sysimage
            ;;
        fallback)
            printf '%s\n' "Warning: sysimage or metadata missing; falling back to plain julia --project=." >&2
            ;;
    esac
}

CURRENT_VERSION=$(get_julia_version)
CURRENT_ARCH=$(get_platform_arch)
CURRENT_HEAD=$(get_head_commit)

if [ ! -f "$SYSIMAGE_PATH" ] || [ ! -f "$META_PATH" ]; then
    handle_missing_sysimage
fi

USE_SYSIMAGE=0
if [ -f "$SYSIMAGE_PATH" ] && [ -f "$META_PATH" ]; then
    META_VERSION=$(json_value julia_version)
    META_FAMILY=$(json_value platform_family)
    META_ARCH=$(json_value platform_arch)
    MISMATCH_REASON=

    if [ -z "$META_VERSION" ]; then
        MISMATCH_REASON="metadata missing julia_version"
    elif [ "$META_VERSION" != "$CURRENT_VERSION" ]; then
        MISMATCH_REASON="Julia version $META_VERSION does not match current Julia $CURRENT_VERSION"
    elif [ -n "$META_FAMILY" ] && [ "$META_FAMILY" != "$PLATFORM_FAMILY" ]; then
        MISMATCH_REASON="platform family $META_FAMILY does not match current platform $PLATFORM_FAMILY"
    elif [ -n "$META_ARCH" ] && [ "$META_ARCH" != "$CURRENT_ARCH" ]; then
        MISMATCH_REASON="platform arch $META_ARCH does not match current arch $CURRENT_ARCH"
    elif [ -n "$CURRENT_HEAD" ] && [ -z "$(json_value git_commit)" ]; then
        MISMATCH_REASON="metadata missing git_commit"
    elif [ -n "$CURRENT_HEAD" ] && [ "$(json_value git_commit)" != "$CURRENT_HEAD" ]; then
        MISMATCH_REASON="sysimage git commit $(json_value git_commit) does not match current HEAD $CURRENT_HEAD"
    fi

    if [ -z "$MISMATCH_REASON" ]; then
        USE_SYSIMAGE=1
    else
        case "$MISMATCH_POLICY" in
            strict)
                printf '%s\n' "Incompatible sysimage: $MISMATCH_REASON" >&2
                exit 1
                ;;
            rebuild)
                printf '%s\n' "Warning: incompatible sysimage detected; rebuilding local sysimage. Reason: $MISMATCH_REASON" >&2
                build_local_sysimage
                META_VERSION=$(json_value julia_version)
                META_FAMILY=$(json_value platform_family)
                META_ARCH=$(json_value platform_arch)
                META_COMMIT=$(json_value git_commit)
                if [ -z "$META_VERSION" ] || [ "$META_VERSION" != "$CURRENT_VERSION" ] || { [ -n "$META_FAMILY" ] && [ "$META_FAMILY" != "$PLATFORM_FAMILY" ]; } || { [ -n "$META_ARCH" ] && [ "$META_ARCH" != "$CURRENT_ARCH" ]; } || { [ -n "$CURRENT_HEAD" ] && { [ -z "$META_COMMIT" ] || [ "$META_COMMIT" != "$CURRENT_HEAD" ]; }; }; then
                    printf '%s\n' "Rebuilt sysimage is still incompatible." >&2
                    exit 1
                fi
                USE_SYSIMAGE=1
                ;;
            fallback)
                printf '%s\n' "Warning: incompatible sysimage detected; falling back to plain julia --project=. Reason: $MISMATCH_REASON" >&2
                ;;
        esac
    fi
fi

if [ "$USE_SYSIMAGE" -eq 1 ]; then
    exec julia --sysimage="$SYSIMAGE_PATH" --project="$REPO_ROOT" "$@"
else
    exec julia --project="$REPO_ROOT" "$@"
fi
