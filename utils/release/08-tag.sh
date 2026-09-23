#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 0 ]] || die "usage: $0"
require_command git
load_release_state
require_marker testpypi-verified
require_clean_worktree
require_artifacts
readonly TAG="v$RELEASE_VERSION"

preflight_remote_artifacts "https://pypi.org/pypi/chargefw/$RELEASE_VERSION/json"
confirm_release_action "Create and push tag $TAG"
if git rev-parse --verify --quiet "refs/tags/$TAG" >/dev/null; then
    [[ $(git cat-file -t "$TAG") == tag ]] || die "local tag $TAG is not annotated"
    [[ $(git rev-list -n 1 "$TAG") == "$RELEASE_COMMIT" ]] ||
        die "local tag $TAG points to another commit"
else
    git tag -a "$TAG" "$RELEASE_COMMIT" -m "ChargeFW $RELEASE_VERSION"
fi
git push origin "refs/tags/$TAG"

remote_matches=false
while IFS=$'\t' read -r object reference; do
    if [[ "$object" == "$RELEASE_COMMIT" && \
          ("$reference" == "refs/tags/$TAG" || "$reference" == "refs/tags/$TAG^{}") ]]; then
        remote_matches=true
    fi
done < <(git ls-remote origin "refs/tags/$TAG" "refs/tags/$TAG^{}")
[[ "$remote_matches" == true ]] || die "remote tag $TAG does not resolve to the release commit"

write_marker tag-pushed
printf 'Created and pushed %s.\n' "$TAG"
