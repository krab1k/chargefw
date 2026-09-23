#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 0 ]] || die "usage: $0"
require_command uv
load_release_state
require_marker testpypi-uploaded
require_artifacts

verify_remote_artifacts "https://test.pypi.org/pypi/chargefw/$RELEASE_VERSION/json"
smoke_test_index "https://test.pypi.org/pypi/chargefw/$RELEASE_VERSION/json" \
    "$RELEASE_ROOT/testpypi-smoke"

write_marker testpypi-verified
printf 'Verified ChargeFW %s from TestPyPI.\n' "$RELEASE_VERSION"
