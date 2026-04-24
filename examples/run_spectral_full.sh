#!/usr/bin/env bash
#==============================================================================
# IOSEE v1.0.0 — run_spectral_full.sh
#
# Full-spectrum run: 10 - 3250 cm^-1 at 1 cm^-1 spacing, linear view-angle
# and wind sweeps (Option 2, via n_view_angle / n_wind).
#
# Intended for benchmarking and for producing a reference grid for
# radiative transfer lookup tables.
#==============================================================================

set -euo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${HERE}"

if   [[ -x "${HERE}/../build/iosee" ]]; then IOSEE="${HERE}/../build/iosee"
elif [[ -x "${HERE}/iosee"          ]]; then IOSEE="${HERE}/iosee"
else
    echo "ERROR: iosee executable not found. Run 'make release' first." >&2
    exit 1
fi

CONFIG="spectral_option2.config"
OUTPUT="emissivity_spectral_full.nc"

echo "=========================================================="
echo "IOSEE full-spectrum spectral run"
echo "  config  : ${CONFIG}    (wavenum 10-3250 cm^-1, 3241 pts)"
echo "  output  : ${OUTPUT}"
echo "=========================================================="

time "${IOSEE}" "${CONFIG}" "output=${OUTPUT}"
