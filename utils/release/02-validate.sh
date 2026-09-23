#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 0 ]] || die "usage: $0"
for command in cmake ctest git pre-commit; do
    require_command "$command"
done
require_clean_worktree

VERSION=$(current_version)
readonly VERSION
validate_version "$VERSION"
if [[ -f "$RELEASE_STATE_FILE" && -f "$(marker_path built)" ]]; then
    # The generated state file is trusted local release state, not user input.
    # shellcheck disable=SC1090
    source "$RELEASE_STATE_FILE"
    if [[ "$RELEASE_VERSION" == "$VERSION" || ! -f "$(marker_path release-complete)" ]]; then
        die "release artifacts already exist; continue that release instead of rebuilding"
    fi
    require_marker release-complete
    readonly PREVIOUS_VERSION=$RELEASE_VERSION
    unset RELEASE_VERSION RELEASE_COMMIT ARTIFACT_DIRECTORY CHECKSUM_FILE
    printf 'Replacing completed release state for ChargeFW %s.\n' "$PREVIOUS_VERSION"
fi
COMMIT=$(git rev-parse HEAD)
readonly COMMIT
initialize_release_state "$VERSION" "$COMMIT"

pre-commit run --all-files

for preset in gcc-debug clang-debug gcc-release clang-release; do
    cmake --preset "$preset"
    cmake --build "build/$preset"
    ctest --test-dir "build/$preset" --output-on-failure -E '^cpptest$'
done

for preset in clang-asan clang-ubsan; do
    cmake --preset "$preset"
    cmake --build "build/$preset" --parallel 1
    ctest --test-dir "build/$preset" --output-on-failure -E '^cpptest$'
done

cmake --preset clang-tidy
cmake --build build/clang-tidy

ENGINE=$(container_engine)
readonly ENGINE
"$ENGINE" build --tag "chargefw:$VERSION" .

require_clean_worktree
[[ $(git rev-parse HEAD) == "$COMMIT" ]] || die "HEAD changed during validation"
load_release_state
write_marker validated
printf 'Validated ChargeFW %s at %s.\n' "$RELEASE_VERSION" "$RELEASE_COMMIT"
