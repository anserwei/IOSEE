#!/usr/bin/env bash
#==============================================================================
# IOSEE v1.0.0 — run_no_multiple_reflection.sh
#
# Demonstrates the 'no_multiple_reflection' flag: outputs the first-order
# Fresnel emissivity (1 - R_F) weighted by the Cox-Munk slope distribution,
# without the multiple-reflection contribution.
#
# This script generates two runs side-by-side (full vs first-order) so the
# user can diff the two NetCDF files with 'ncdump' or 'cdo diff'.
#
# Note: the shipped config files do not set no_multiple_reflection. This
# script patches a working copy on the fly (leaving the originals intact).
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

BASE_CONFIG="spectral_option1.config"
FULL_CONFIG="_run_full.config"
NOMR_CONFIG="_run_nomr.config"

trap 'rm -f "${FULL_CONFIG}" "${NOMR_CONFIG}"' EXIT

# full (default) run - copy config verbatim
cp "${BASE_CONFIG}" "${FULL_CONFIG}"

# first-order-only run - inject no_multiple_reflection = .true.
# inserted right after the opening &configuration line
awk 'NR==1{print; print "no_multiple_reflection = .true.,"; next} 1' \
    "${BASE_CONFIG}" > "${NOMR_CONFIG}"

echo "=========================================================="
echo "IOSEE: full emissivity (with multiple reflection)"
echo "=========================================================="
"${IOSEE}" "${FULL_CONFIG}" "output=emiss_full.nc"

echo ""
echo "=========================================================="
echo "IOSEE: first-order emissivity (no multiple reflection)"
echo "=========================================================="
"${IOSEE}" "${NOMR_CONFIG}" "output=emiss_nomr.nc"
