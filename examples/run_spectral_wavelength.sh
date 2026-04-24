#!/usr/bin/env bash
#==============================================================================
# IOSEE v1.0.0 — run_spectral_wavelength.sh
#
# Wavelength input mode example. The driver automatically detects the
# wavelen_* flags in the namelist and converts internally to wavenumber.
#
# Config:  spectral_option1_wl.config  (8 - 14 um window used)
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

CONFIG="spectral_option1_wl.config"
OUTPUT="emissivity_spectral_wavelength.nc"

echo "=========================================================="
echo "IOSEE wavelength-input spectral run"
echo "  config  : ${CONFIG}"
echo "  output  : ${OUTPUT}"
echo "=========================================================="

"${IOSEE}" "${CONFIG}" "output=${OUTPUT}"
