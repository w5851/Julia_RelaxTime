#!/usr/bin/env sh

set -eu

REPO="w5851/Julia_RelaxTime"
RELEASE_TAG="latest"
OVERWRITE=0
DRY_RUN=0

while [ "$#" -gt 0 ]; do
    case "$1" in
        --repo=*)
            REPO=${1#*=}
            shift
            ;;
        --release-tag=*)
            RELEASE_TAG=${1#*=}
            shift
            ;;
        --overwrite)
            OVERWRITE=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --help|-h)
            cat <<'EOF'
Usage:
  scripts/dev/bootstrap_sysimage.sh [--repo=<owner/name>] [--release-tag=<tag|latest>] [--overwrite] [--dry-run]

Behavior:
  - resolve current platform family, arch, and Julia version
  - construct the matching release asset name
  - download and extract the prebuilt sysimage into build/
EOF
            exit 0
            ;;
        *)
            printf '%s\n' "Unsupported argument: $1" >&2
            exit 1
            ;;
    esac
done

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
REPO_ROOT=$(CDPATH= cd -- "$SCRIPT_DIR/../.." && pwd)
OUTPUT_DIR="$REPO_ROOT/build"

uname_s=$(uname -s 2>/dev/null || printf '')
uname_m=$(uname -m 2>/dev/null || printf '')

case "$uname_s" in
    Darwin) PLATFORM_FAMILY=macos ; SYSIMAGE_EXT=dylib ; ARCHIVE_FORMAT=tar.gz ;;
    Linux) PLATFORM_FAMILY=linux ; SYSIMAGE_EXT=so ; ARCHIVE_FORMAT=tar.gz ;;
    MINGW*|MSYS*|CYGWIN*) PLATFORM_FAMILY=windows ; SYSIMAGE_EXT=dll ; ARCHIVE_FORMAT=zip ;;
    *)
        printf '%s\n' "Unsupported platform for sysimage bootstrap: $uname_s" >&2
        exit 1
        ;;
esac

case "$uname_m" in
    amd64|x64) PLATFORM_ARCH=x86_64 ;;
    arm64) PLATFORM_ARCH=aarch64 ;;
    *) PLATFORM_ARCH=$uname_m ;;
esac

JULIA_VERSION=$(julia --version | sed -E 's/.* ([0-9]+\.[0-9]+\.[0-9]+).*/\1/')
ASSET_NAME="JuliaRelaxTime-sysimage-$PLATFORM_FAMILY-$PLATFORM_ARCH-julia$JULIA_VERSION.$ARCHIVE_FORMAT"

if [ "$RELEASE_TAG" = "latest" ]; then
    ASSET_URL="https://github.com/$REPO/releases/latest/download/$ASSET_NAME"
else
    ASSET_URL="https://github.com/$REPO/releases/download/$RELEASE_TAG/$ASSET_NAME"
fi

SYSIMAGE_TARGET="$OUTPUT_DIR/JuliaRelaxTime.$SYSIMAGE_EXT"
META_TARGET="$OUTPUT_DIR/JuliaRelaxTime.sysimage.json"

if [ "$OVERWRITE" -ne 1 ] && { [ -f "$SYSIMAGE_TARGET" ] || [ -f "$META_TARGET" ]; }; then
    printf '%s\n' "Target sysimage or metadata already exists in $OUTPUT_DIR. Re-run with --overwrite to replace it." >&2
    exit 1
fi

if [ "$DRY_RUN" -eq 1 ]; then
    printf 'Resolved asset: %s\n' "$ASSET_NAME"
    printf 'Download URL: %s\n' "$ASSET_URL"
    printf 'Output dir: %s\n' "$OUTPUT_DIR"
    exit 0
fi

mkdir -p "$OUTPUT_DIR"
TMP_ROOT=$(mktemp -d "${TMPDIR:-/tmp}/jrt-sysimage.XXXXXX")
ARCHIVE_PATH="$TMP_ROOT/$ASSET_NAME"
EXTRACT_DIR="$TMP_ROOT/extract"
mkdir -p "$EXTRACT_DIR"

cleanup() {
    rm -rf "$TMP_ROOT"
}
trap cleanup EXIT INT TERM

if command -v curl >/dev/null 2>&1; then
    curl -fL "$ASSET_URL" -o "$ARCHIVE_PATH"
elif command -v wget >/dev/null 2>&1; then
    wget -O "$ARCHIVE_PATH" "$ASSET_URL"
else
    printf '%s\n' "Neither curl nor wget is available for download." >&2
    exit 1
fi

if [ "$ARCHIVE_FORMAT" = "zip" ]; then
    if command -v unzip >/dev/null 2>&1; then
        unzip -q "$ARCHIVE_PATH" -d "$EXTRACT_DIR"
    else
        printf '%s\n' "unzip is required to extract zip assets." >&2
        exit 1
    fi
else
    tar -xzf "$ARCHIVE_PATH" -C "$EXTRACT_DIR"
fi

SYSIMAGE_FILE=$(find "$EXTRACT_DIR" -type f -name "JuliaRelaxTime.$SYSIMAGE_EXT" | head -n 1)
META_FILE=$(find "$EXTRACT_DIR" -type f -name "JuliaRelaxTime.sysimage.json" | head -n 1)

if [ -z "$SYSIMAGE_FILE" ] || [ -z "$META_FILE" ]; then
    printf '%s\n' "Release asset does not contain JuliaRelaxTime.$SYSIMAGE_EXT and JuliaRelaxTime.sysimage.json" >&2
    exit 1
fi

cp "$SYSIMAGE_FILE" "$SYSIMAGE_TARGET"
cp "$META_FILE" "$META_TARGET"

printf 'Bootstrapped sysimage: %s\n' "$SYSIMAGE_TARGET"
printf 'Bootstrapped metadata: %s\n' "$META_TARGET"
