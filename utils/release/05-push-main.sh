#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 0 ]] || die "usage: $0"
require_command git
load_release_state
require_marker local-verified
require_clean_worktree
[[ $(git branch --show-current) == main ]] ||
    die "switch or fast-forward main to $RELEASE_COMMIT before running this step"

confirm_release_action "Push the release commit to origin/main"
git push origin main
REMOTE_COMMIT=$(git ls-remote origin refs/heads/main | cut -f1)
readonly REMOTE_COMMIT
[[ "$REMOTE_COMMIT" == "$RELEASE_COMMIT" ]] || die "origin/main does not point to the release commit"

write_marker main-pushed
printf 'Pushed ChargeFW %s commit to origin/main.\n' "$RELEASE_VERSION"
