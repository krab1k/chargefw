#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 0 ]] || die "usage: $0"
for command in git uv; do
    require_command "$command"
done
load_release_state
require_marker validated
require_clean_worktree
[[ ! -f "$(marker_path built)" ]] ||
    die "release artifacts are already built; do not rebuild them"

rm -rf -- "$ARTIFACT_DIRECTORY"
mkdir -p -- "$ARTIFACT_DIRECTORY"
uv build --sdist --out-dir "$ARTIFACT_DIRECTORY" --clear

shopt -s nullglob
sdists=("$ARTIFACT_DIRECTORY"/chargefw-*.tar.gz)
[[ ${#sdists[@]} -eq 1 ]] || die "expected exactly one source distribution"
[[ $(basename -- "${sdists[0]}") == "chargefw-$RELEASE_VERSION.tar.gz" ]] ||
    die "source distribution has an unexpected version"

ENGINE=$(container_engine)
readonly ENGINE
CIBW_CONTAINER_ENGINE="$ENGINE" uvx cibuildwheel==3.4.1 \
    --platform linux \
    --output-dir "$ARTIFACT_DIRECTORY" \
    "${sdists[0]}"

{
    printf 'version: %s\n' "$RELEASE_VERSION"
    printf 'commit: %s\n' "$RELEASE_COMMIT"
    printf 'built-at-utc: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    git --version
    uv --version
    "$ENGINE" --version
} > "$RELEASE_ROOT/build-info.txt"

write_marker built
printf 'Built release artifacts in %s.\n' "$ARTIFACT_DIRECTORY"
