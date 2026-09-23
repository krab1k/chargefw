#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 0 ]] || die "usage: $0"
require_command uv
load_release_state
require_marker main-pushed
require_artifacts
[[ ${TWINE_USERNAME:-} == __token__ ]] || die "set TWINE_USERNAME=__token__"
[[ -n ${TWINE_PASSWORD:-} ]] || die "set TWINE_PASSWORD to a TestPyPI project token"

confirm_release_action "Upload release artifacts to TestPyPI"
upload_missing_artifacts "https://test.pypi.org/legacy/" \
    "https://test.pypi.org/pypi/chargefw/$RELEASE_VERSION/json"

write_marker testpypi-uploaded
printf 'Uploaded ChargeFW %s artifacts to TestPyPI.\n' "$RELEASE_VERSION"
