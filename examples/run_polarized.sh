#!/usr/bin/env bash
#==============================================================================
# IOSEE v1.0.0 — run_polarized.sh
#
# V- and H-polarized emissivity (Henderson et al. 2003) in the
# 8 - 12 um atmospheric window.
#
# Config:  spectral_option1_polarized_wl.config
# Output:  NetCDF with emissivity, emissivity_v, emissivity_h.
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

CONFIG="spectral_option1_polarized_wl.config"
OUTPUT="emissivity_polarized_8-12um.nc"

echo "=========================================================="
echo "IOSEE polarized spectral run (V and H)"
echo "  config  : ${CONFIG}"
echo "  output  : ${OUTPUT}"
echo "=========================================================="

"${IOSEE}" "${CONFIG}" "output=${OUTPUT}"
