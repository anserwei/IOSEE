# IOSEE Changelog

All notable capability changes to the IOSEE (Infrared Ocean Surface
Effective Emissivity) package are documented in this file.

Only capability additions, improvements, and interface changes are
recorded. Internal maintenance and routine bug fixes are not included.

Format based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and [Semantic Versioning](https://semver.org/).

> **Note on version numbering.** Versions `0.8.1` through `0.9.6`
> constituted the public beta of IOSEE. Release **1.0.0** is the first
> stable release with a frozen public API and the full physics pipeline
> validated across the 10 – 5000 cm⁻¹ infrared range. The previous
> development tags `1.0.0 … 1.0.12` have been renumbered to the
> semver-compliant `0.8.1 … 1.0.0` scheme shown below.

---

## [1.0.0] — 2026-03-21 — Stable Release

### Summary

First stable release. Consolidates the 0.8.1 → 0.9.6 beta series and
freezes the public API. Comment and documentation cleanup, updated
Wang et al. (2026) citation (published), and replaced "observational"
with neutral terminology in code comments and physics documentation.

### Changed

- `lib/mathlib.F90` — `smooth_effective_angle_lut()` comments:
  "observational anchors" → "reference anchors".
- `src/sea_water_nk.F90` — citation:
  "Wang et al. (2025)" → "Wang et al. (2026)" in module header and
  inline comments.
- `docs/physics/physics_deep_dive.md` — terminology and citation:
  "observational anchors" → "reference anchors";
  "tuning algorithm" → "merging algorithm";
  full Wang et al. (2026) reference added.
- `CHANGELOG.txt` → `CHANGELOG.md`; `VERSION.txt` → `VERSION.md`
  (Markdown format).

### Documentation

- `docs/user_guide/configurable_flag.docx` — new, comprehensive table
  of every namelist flag with valid options and defaults.
- `examples/*.sh` — new shell wrappers demonstrating canonical
  run patterns (spectral, narrow-band, polarized, wavelength-input,
  parallel-thread control).

### Data Source

- `data/water_optical_constants.nc`:

  > Wang, S., Yang, P., Brindley, H. E., Huang, X., & L'Ecuyer, T. S.
  > (2026). *Enhanced full spectral temperature-dependent refractive
  > index of liquid water from supercooled to ambient conditions.*
  > Geophysical Research Letters, 53(6), e2025GL119385.
  > https://doi.org/10.1029/2025GL119385

---

## [0.9.6] — 2026-03-04

### Summary

Documentation improvements and codebase housekeeping.

### Documentation

- `docs/user_guide/configuration.md`:
  - Comprehensive documentation of all configuration parameters.
  - Both wavenumber and wavelength input methods described.
  - Explicit list and range specification for view angles and winds.
  - All flags with valid ranges, defaults, and physical meaning.
  - Annotated examples for each calculation mode.
- `NOTICE` — added attribution for ECMWF-derived NetCDF I/O utilities.

---

## [0.9.5] — 2026-02-26

### Summary

Wind-adaptive blending for the effective view angle LUT smoothing. At
high wind speeds (> 10 m/s) the effective-to-geometric angle gap grows
dramatically (≈ 24° at 14 m/s), producing a visible emissivity kink at
70°. The blend window now adapts from [60°, 70°] at calm winds to
[50°, 70°] at 20 m/s, spreading the transition and eliminating the kink.

### Changed

- `lib/mathlib.F90` — `smooth_effective_angle_lut()`:
  - Subroutine now accepts `wind_table(:)` for wind-adaptive blending.
  - Phase-2 `blend_start` adapts per wind speed via cosine interpolation:

    ```
    wind ≤ 10 m/s:  blend_start = 60°  (unchanged from 0.9.4)
    wind = 20 m/s:  blend_start = 50°  (20° window)
    10 < wind < 20: blend_start = 60 − w_wind · 10
                    w_wind = 0.5·(1 − cos(π·(wind − 10)/10))
    ```

  - C¹ continuity preserved at both endpoints for all wind speeds.
  - Wind speeds ≤ 10 m/s produce bit-for-bit identical output to 0.9.4.

- `src/driver.F90` — effective-view-angle loading section:
  passes `reference_wind_speed` to `smooth_effective_angle_lut()`.

### Documentation

- `docs/physics/physics_deep_dive.md` — updated Phase-2 section with the
  wind-adaptive blend formula and per-wind constraint preservation table.

### Physical Background

At high wind speeds the effective-to-geometric gap widens:

| Angle | Eff(wind = 0) | Eff(wind = 14) | Gap   |
|-------|---------------|----------------|-------|
| 55°   | 51.7°         | 38.1°          | 16.9° |
| 60°   | 58.1°         | 46.2°          | 13.8° |
| 65°   | 62.9°         | 55.1°          |  9.9° |
| 70°   | 63.5°         | 63.0°          |  7.0° |

The wind-adaptive window spreads the transition to ≈ 13.5° at 14 m/s,
reducing the peak effective-angle rate from ≈ 2.4 deg/deg to ≈ 1.8 deg/deg.

### Backward Compatibility

- Wind speeds 0–10 m/s: bit-for-bit identical to 0.9.4.
- Wind speeds > 10 m/s: smoother transitions at the 70° boundary.
- No API changes; configuration files work without modification.
- Output format unchanged (same NetCDF structure).

---

## [0.9.4] — 2026-02-25

### Summary

Smooth emissivity curves at the 70° effective/geometric angle boundary.
Eliminates the emissivity discontinuity at 70° caused by the effective
view angle LUT plateau, using PCHIP smoothing of the effective angle LUT.

### Added

#### Effective view angle LUT smoothing

- `lib/mathlib.F90` — `smooth_effective_angle_lut()` subroutine:

  Three-phase smoothing applied to the effective view angle LUT at
  runtime:

  1. **Phase 1** — PCHIP re-interpolation through 5° reference anchors.
  2. **Phase 2** — blend toward the geometric angle over [60°, 70°].
  3. **Phase 3** — non-negativity and monotonicity enforcement.

  Applied to all wind speeds; wavenumber-independent (LUT is
  angle × wind).

- `src/driver.F90` — `enforce_emissivity_monotonicity()` subroutine:
  forward-clamping safety net for any residual non-monotonicity in the
  final emissivity-vs-angle curves.

### Changed

- `src/driver.F90` — effective-view-angle loading section:
  calls `smooth_effective_angle_lut()` after loading the LUT from NetCDF.

### Documentation

- `docs/physics/physics_deep_dive.md` — new Section 6 "Effective View
  Angle Processing" covering the full pipeline:
  - LUT generation via spectral variance minimization
    (Nalli et al. 2008 / 2023).
  - Three-phase LUT smoothing algorithm with mathematical details.
  - 70° boundary blending and monotonicity enforcement.
  - End-to-end data-flow diagram.
- `docs/physics/README.md` — updated to reference Section 6.
- `docs/architecture/architecture_overview.md` — updated mathlib section.

### Physical Background

The effective view angle LUT (Nalli et al. 2008, 2023) maps geometric
viewing angles (0–70°) to effective incidence angles that account for
wave-slope statistics. At 70° the raw LUT shows a plateau
(eff_angle ≈ 63° for all wind speeds 0–20 m/s), a boundary artifact of
the spectral-variance-minimization fit.

When the driver transitions from effective to geometric angle at 70°,
the Fresnel emissivity jumps from F(63°) to F(70°) — a 37× change in
the emissivity derivative. The smoothing forces the effective angle to
equal the geometric angle at exactly 70°, eliminating the discontinuity.

### Backward Compatibility

- User angles ≤ 70°: minor emissivity changes in 60–70° from smoothing.
- User angles > 70°: the 70° emissivity discontinuity is eliminated.
- No API changes; configuration files work without modification.
- Output format unchanged (same NetCDF structure).

---

## [0.9.3] — 2026-02-25

### Summary

Extended viewing angle support to 85° via angle blending at the 70°
boundary. Angles ≤ 70° use effective view angles (Nalli et al.);
angles > 70° use normal (geometric) view angles with a continuity
correction at the 70° boundary.

### Added

#### Angle blending at the 70° boundary

- `src/driver.F90` — automatic angle blending for viewing angles > 70°:
  - Loads `effective_view_angle.nc` (0–70°) and
    `normal_view_angle.nc` (0–85°).
  - Splits user angles at 70° into two groups:
    * Angles ≤ 70° → emissivity via effective view angles.
    * Angles > 70° → emissivity via normal (geometric) view angles.
  - Applies continuity correction at boundary:
    ```
    diff_70 = emiss_eff(70°) − emiss_nor(70°)
    emiss_final(>70°) = emiss_nor(>70°) + diff_70
    ```
  - Adds 70° as bridge point automatically if absent from user's angles.
  - Falls back to original behavior when all angles ≤ 70° (no overhead).

- `src/driver.F90` — `blend_angle_results()` helper subroutine:
  generic blending routine used by all processing modes; computes
  per-wind, per-wavenumber continuity correction and maps effective and
  normal partial results into the final output array.

#### Supported in all processing modes

- Spectral mode (unpolarized), including `no_multiple_reflection`.
- Polarized spectral mode: blends V, H, no-reflection, and reflectance.
- Narrow-band mode: blending applied on the full spectrum before
  Planck-function integration.

### Physical Background

The effective view angle correction (Nalli et al. 2008, 2023) accounts
for the difference between the satellite viewing geometry and the
effective incidence angle due to wave-slope statistics at the rough
ocean surface. This correction is valid up to 70°. Beyond 70°, the
physical (geometric) viewing angle is used. The continuity correction
ensures smooth emissivity transitions across the boundary:

```
emiss_eff(0–70°)  : computed with effective view angles from LUT
emiss_nor(70–85°) : computed with geometric (normal) view angles
diff_70 = emiss_eff(70°) − emiss_nor(70°)
emiss_final(≤70°) = emiss_eff
emiss_final(>70°) = emiss_nor + diff_70
```

### Backward Compatibility

- All user angles ≤ 70°: behavior is **identical** to 0.9.2.
- Angles > 70° present: new blending logic activates automatically.
- No API changes; configuration files work without modification.
- Output format unchanged (same NetCDF structure).

---

## [0.9.2] — 2026-01-25

### Summary

Removed the legacy `get_emissivity` subroutine. The public API now
exposes only the optimized interfaces introduced in 0.9.0.

### Removed

- `src/ocean_emissivity.F90` — `get_emissivity` subroutine (≈ 127 lines):
  superseded by `get_emissivity_optimized`. Removal simplifies the
  public API and reduces maintenance burden.

### Changed

Public interface in `ocean_emissivity.F90` now consists of:

- `get_emissivity_optimized` — optimized unified wind processing
  (**primary** entry point).
- `get_polarized_emissivity_optimized` — V/H polarized calculations.
- `initialize_workspace_pool_emiss` — thread-safe workspace init.
- `validate_input_dimensions` — input validation utility.

---

## [0.9.1] — 2026-01-24

### Summary

New configuration option `no_multiple_reflection` for output control.
When enabled, outputs first-order Fresnel emissivity without
multiple-reflection contributions (`emiss_norefle`) instead of the full
emissivity. An internal Fresnel LUT range extension (`n_max` 1.6 → 3.0,
`k_max` 0.5 → 1.5) was also folded into this release.

### Added

#### New configuration option

- `lib/configuration.F90` — `no_multiple_reflection`:
  - New `logical` parameter (default: `.false.`).
  - When `.true.`, outputs emissivity **without** multiple reflection
    effects.
  - Physical meaning: first-order emissivity
    `(1 − Fresnel_reflectance)` weighted over the Cox-Munk slope
    distribution (Masuda 2006).
  - Available in both spectral and narrow-band modes.

#### Fresnel LUT range extension

- `src/ocean_emissivity.F90` — Fresnel LUT extended:
  `n_max` 1.6 → 3.0, `k_max` 0.5 → 1.5 for robust coverage of the full
  10 – 5000 cm⁻¹ range.

#### Config file usage

```fortran
&configuration
no_multiple_reflection = .true.,
/
```

### Output Behavior

| `no_multiple_reflection` | NetCDF `emissivity` variable contains            |
|--------------------------|--------------------------------------------------|
| `.false.` (default)      | Full emissivity with multiple reflection         |
| `.true.`                 | First-order emissivity (no multiple reflection)  |

### Physical Background

The ocean surface emissivity is

```
emiss_final = emiss_norefle + multi_refle
```

where `emiss_norefle` is the first-order emissivity
(`1 − Fresnel_reflectance`) weighted over the Cox-Munk slope
distribution and `multi_refle` is the additional contribution from
multiple reflections of atmospheric downwelling radiation between wave
facets.

### Backward Compatibility

- Default behavior unchanged (`no_multiple_reflection = .false.`).
- Existing configuration files work without modification.
- Output format unchanged.

---

## [0.9.0] — 2026-01-24

### Summary

Performance optimization: unified wind-processing architecture.
Eliminates the serial wind loop bottleneck by matching the
high-performance architecture of `get_polarized_emissivity_optimized`.
All computational algorithms remain identical to 0.8.3.

### Added

#### New optimized subroutine

- `src/ocean_emissivity.F90` — `get_emissivity_optimized`:
  - Accepts a 2-D `effective_angles_2d(:,:)` array directly (no
    extraction overhead).
  - Pre-allocated output arrays with `intent(inout)`.
  - Single function call processes **all** wind conditions at once.
  - `COLLAPSE(3)` parallelization across all dimensions
    (angles × winds × wavenumbers).
  - `SCHEDULE(guided)` for optimal dynamic load balancing.
  - Core algorithm identical to the original `get_emissivity`.

### Changed

- `src/driver.F90` — `process_spectral_mode()`: serial wind loop
  replaced with a single `get_emissivity_optimized` call.
- `src/driver.F90` — `process_narrow_band_mode()`: serial wind loop
  replaced with a single `get_emissivity_optimized` call; pre-allocates
  `full_emiss(n_angles, n_winds, full_n_points)` for all winds.

### Performance

| Mode        | Winds | Previous     | Optimized      | Speedup |
|-------------|-------|--------------|----------------|---------|
| Spectral    | 5     | 12 500 c/s   | ≈ 18 000 c/s   | ≈ 1.4×  |
| Spectral    | 10    | 12 500 c/s   | ≈ 22 000 c/s   | ≈ 1.8×  |
| Narrow-band | any   | variable     | consistent     | 2–4×    |

### Backward Compatibility

- Public API unchanged for external callers.
- Output format unchanged (NetCDF files are identical).
- Numerical results are bit-for-bit identical to 0.8.3.

---

## [0.8.3] — 2026-01-24

### Summary

Multi-compiler and multi-platform support with bit-for-bit
reproducibility option. Adds Intel ifx (LLVM) and Cray ftn compilers and
automatic parallel build detection. All computational algorithms remain
identical to 0.8.2.

### Added

#### New compiler support

- `cmake/CompilerFlagsIntelLLVM.cmake`:
  - Full support for Intel ifx (2024.0+) with LLVM backend.
  - Platform-aware optimization: `-xHost` on workstations,
    `-xCORE-AVX2` on HPC.
  - Strict IEEE 754 mode with
    `-fp-model=strict -fimf-arch-consistency=true`.

- `cmake/CompilerFlagsCray.cmake`:
  - Full support for Cray Fortran (ftn) from CCE 15.0+.
  - Uses Cray-specific flags: `-h omp` for OpenMP, `-h fp0` for
    reproducibility.
  - Automatic detection of `CRAY_CPU_TARGET` environment variable.

#### Reproducibility mode

- CMake option `ENABLE_REPRODUCIBLE` (default: `OFF`):
  - Produces bit-for-bit identical results across different compilers
    and platforms.
  - Strict IEEE 754 compliance with no fast-math optimizations.
  - Build with: `cmake -DENABLE_REPRODUCIBLE=ON ..`
    or `make reproducible`.

#### Parallel build automation

- Top-level `Makefile`:
  - Auto-detects CPU cores for parallel compilation.
  - Ninja generator support for faster builds.
  - Targets: `make release`, `make debug`, `make reproducible`,
    `make serial`, `make info`.

- `run_iosee.sh` — runtime optimization wrapper:
  - Automatically configures OpenMP environment variables.
  - Detects physical cores (not hyperthreads) for optimal thread count.
  - NUMA-aware execution on multi-socket Linux systems.

### Compiler Support Matrix

| Compiler     | CMake ID   | Min. Version | Reproducibility Flag              |
|--------------|------------|--------------|-----------------------------------|
| GNU gfortran | GNU        | 12.0         | `-ffp-contract=off -fno-fast-math`|
| Intel ifort  | Intel      | 19.0         | `-fp-model strict`                |
| Intel ifx    | IntelLLVM  | 2024.0       | `-fp-model=strict`                |
| Cray ftn     | Cray       | CCE 15       | `-h fp0`                          |

### Backward Compatibility

- All algorithms unchanged; numerical results identical to 0.8.2.
- Configuration files work without modification.

---

## [0.8.2] — 2025-10-24

### Summary

Package renamed from **OCIOSE** to **IOSEE** (Infrared Ocean Surface
Effective Emissivity).

### Changed

- **Renamed** package from `OCIOSE` to `IOSEE`.
- **Renamed** executable from `ocean_emissivity` to `iosee`:
  - Usage: `./iosee config_file.config`.

---

## [0.8.1] — 2025-10-09 — Initial Beta

### Summary

Initial beta release of IOSEE (Infrared Ocean Surface Effective
Emissivity), under the prior project name OCIOSE.

### Features

- Ocean surface emissivity calculations for the infrared spectrum
  (10 – 5000 cm⁻¹).
- Three operational modes:
  - **Spectral** — high-resolution wavenumber calculations.
  - **Narrow-band** — band-averaged emissivity with Planck-function
    weighting.
  - **Polarized** — V/H polarization emissivity output.
- Physics implementations based on peer-reviewed literature:
  - Wu & Smith (1997): emissivity of the rough sea surface.
  - Masuda (2006): multiple-reflection effects with Cox-Munk slope
    distribution.
  - Henderson et al. (2003): polarized emissivity modeling.
  - Nalli et al. (2008, 2023): effective view angle corrections.
- Complex refractive index with temperature and salinity dependence
  (Wang et al. 2026; Röttgers et al. 2014).
- OpenMP parallelization for multi-core performance.
- NetCDF input/output.
- Configurable via Fortran namelist configuration files.

### Dependencies

- CMake ≥ 3.16
- Fortran 2008 compiler (gfortran 12+ or Intel ifx 2024+).
- NetCDF-Fortran library.
- OpenMP (optional).

### Authors

- Dr. Jian Wei (anser@tamu.edu)
- Atmospheric & Oceanic Optics Group, Department of Atmospheric
  Sciences, Texas A&M University.

### License

Apache License, Version 2.0.
