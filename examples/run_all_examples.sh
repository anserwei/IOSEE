#!/usr/bin/env bash
#==============================================================================
# IOSEE v1.0.0 — run_all_examples.sh
#
# Convenience driver: runs every shipped configuration case in sequence.
# Each case writes to a dedicated output NetCDF so results are not
# overwritten. Failed cases print an error banner and the script continues
# to the next.
#==============================================================================

set -u

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "${HERE}"

if   [[ -x "${HERE}/../build/iosee" ]]; then IOSEE="${HERE}/../build/iosee"
elif [[ -x "${HERE}/iosee"          ]]; then IOSEE="${HERE}/iosee"
else
    echo "ERROR: iosee executable not found. Run 'make release' first." >&2
    exit 1
fi

CONFIGS=(
    "spectral_option1.config"
    "spectral_option2.config"
    "spectral_option1_wl.config"
    "spectral_option2_wl.config"
    "spectral_option1_polarized_wn.config"
    "spectral_option1_polarized_wl.config"
    "spectral_option2_polarized_wn.config"
    "spectral_option2_polarized_wl.config"
    "spectral_option2_rrtmg.config"
    "narrow_band_option1.config"
    "narrow_band_option2.config"
    "narrow_band_option2_rrtmg.config"
)

pass=0
fail=0
for cfg in "${CONFIGS[@]}"; do
    base="${cfg%.config}"
    out="emiss_${base}.nc"
    echo "------------------------------------------------------------"
    echo "Running ${cfg} -> ${out}"
    echo "------------------------------------------------------------"
    if "${IOSEE}" "${cfg}" "output=${out}"; then
        pass=$((pass + 1))
    else
        fail=$((fail + 1))
        echo "!!! FAILED: ${cfg}" >&2
    fi
done

echo ""
echo "============================================================"
echo "Completed: ${pass} passed, ${fail} failed out of ${#CONFIGS[@]}"
echo "============================================================"
exit $(( fail > 0 ? 1 : 0 ))
