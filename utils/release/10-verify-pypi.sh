#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 0 ]] || die "usage: $0"
require_command uv
load_release_state
require_marker pypi-uploaded
require_artifacts

verify_remote_artifacts "https://pypi.org/pypi/chargefw/$RELEASE_VERSION/json"
smoke_test_index "https://pypi.org/pypi/chargefw/$RELEASE_VERSION/json" \
    "$RELEASE_ROOT/pypi-smoke"

write_marker release-complete
printf 'ChargeFW %s is published and verified on PyPI.\n' "$RELEASE_VERSION"
