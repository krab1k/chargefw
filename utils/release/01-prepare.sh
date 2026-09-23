#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 1 ]] || die "usage: $0 VERSION"
readonly NEW_VERSION=$1
validate_version "$NEW_VERSION"
require_command git
require_command perl
require_clean_worktree

OLD_VERSION=$(current_version)
readonly OLD_VERSION
[[ -n "$OLD_VERSION" ]] || die "could not read the current version from CMakeLists.txt"
[[ "$NEW_VERSION" != "$OLD_VERSION" ]] || die "CMakeLists.txt already contains version $NEW_VERSION"
if git rev-parse --verify --quiet "refs/tags/v$NEW_VERSION" >/dev/null; then
    die "tag v$NEW_VERSION already exists"
fi

OLD_VERSION="$OLD_VERSION" NEW_VERSION="$NEW_VERSION" perl -0pi -e '
    my $old = quotemeta($ENV{OLD_VERSION});
    my $new = $ENV{NEW_VERSION};
    my $count = s/(project\s*\(\s*chargefw\s+VERSION\s+)$old\b/${1}$new/gis;
    die "expected exactly one project version\n" unless $count == 1;
' CMakeLists.txt

[[ $(current_version) == "$NEW_VERSION" ]] || die "version update failed"
printf 'Updated ChargeFW %s -> %s. Review and commit CMakeLists.txt, then run 02-validate.sh.\n' \
    "$OLD_VERSION" "$NEW_VERSION"
