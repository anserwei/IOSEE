#!/usr/bin/env bash
#==============================================================================
# IOSEE v1.0.0 — run_narrow_band.sh
#
# Narrow-band (Planck-weighted band-averaged) emissivity using the
# RRTMGP longwave band boundaries.
#
# Config:  narrow_band_option2_rrtmg.config
# Output:  one band-averaged emissivity per RRTMGP band.
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

CONFIG="narrow_band_option2_rrtmg.config"
OUTPUT="emissivity_rrtmg_bands.nc"

echo "=========================================================="
echo "IOSEE narrow-band (RRTMGP LW bands) run"
echo "  config  : ${CONFIG}"
echo "  output  : ${OUTPUT}"
echo "=========================================================="

"${IOSEE}" "${CONFIG}" "output=${OUTPUT}"
