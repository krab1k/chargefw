#!/usr/bin/env bash

set -euo pipefail
# shellcheck source=utils/release/_common.sh
source "$(dirname -- "${BASH_SOURCE[0]}")/_common.sh"

[[ $# -eq 0 ]] || die "usage: $0"
require_command uv
load_release_state
require_marker tag-pushed
require_artifacts
[[ ${TWINE_USERNAME:-} == __token__ ]] || die "set TWINE_USERNAME=__token__"
[[ -n ${TWINE_PASSWORD:-} ]] || die "set TWINE_PASSWORD to a PyPI project token"

confirm_release_action "Upload the TestPyPI-approved artifacts to production PyPI"
upload_missing_artifacts "https://upload.pypi.org/legacy/" \
    "https://pypi.org/pypi/chargefw/$RELEASE_VERSION/json"

write_marker pypi-uploaded
printf 'Uploaded ChargeFW %s artifacts to production PyPI.\n' "$RELEASE_VERSION"
