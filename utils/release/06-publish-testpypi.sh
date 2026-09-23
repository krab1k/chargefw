#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 0 ]] || die "usage: $0"
require_command uv
load_release_state
require_marker main-pushed
require_artifacts
[[ ${TWINE_USERNAME:-__token__} == __token__ ]] || die "TWINE_USERNAME must be __token__"

confirm_release_action "Upload release artifacts to TestPyPI"
upload_missing_artifacts testpypi "https://test.pypi.org/legacy/" \
    "https://test.pypi.org/pypi/chargefw/$RELEASE_VERSION/json"

write_marker testpypi-uploaded
printf 'Uploaded ChargeFW %s artifacts to TestPyPI.\n' "$RELEASE_VERSION"
