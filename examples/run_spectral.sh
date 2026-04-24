#!/usr/bin/env bash
#==============================================================================
# IOSEE v1.0.0 — run_spectral.sh
#
# Canonical example: high-resolution spectral emissivity on a single band.
#
# Uses:
#   mode       = 'spectral'
#   wavenum    = 530 - 1250 cm^-1, 145 points
#   view_angle = 25..75 deg (explicit list, Option 1)
#   wind       = 0..18 m/s  (explicit list, Option 1)
#
# Run from the top-level repository directory:
#     bash examples/run_spectral.sh
#
# Or from inside ./examples (uses the iosee symlink produced by cmake):
#     cd examples && bash run_spectral.sh
#==============================================================================

set -euo pipefail

# Locate the iosee executable (prefer ../build/iosee, fall back to ./iosee)
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${HERE}"

if   [[ -x "${HERE}/../build/iosee" ]]; then IOSEE="${HERE}/../build/iosee"
elif [[ -x "${HERE}/iosee"          ]]; then IOSEE="${HERE}/iosee"
else
    echo "ERROR: iosee executable not found. Run 'make release' first." >&2
    exit 1
fi

CONFIG="spectral_option1.config"
OUTPUT="emissivity_spectral_option1.nc"

echo "=========================================================="
echo "IOSEE spectral run"
echo "  exe     : ${IOSEE}"
echo "  config  : ${CONFIG}"
echo "  output  : ${OUTPUT}"
echo "=========================================================="

"${IOSEE}" "${CONFIG}" "output=${OUTPUT}"
