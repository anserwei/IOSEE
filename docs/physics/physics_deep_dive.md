# IOSEE Physics Deep-Dive

This document provides a detailed technical analysis of the physical models and their implementation in IOSEE.

---

## 1. Fresnel Equations Implementation

### 1.1 Physical Foundation

The Fresnel equations describe the reflection and transmission of electromagnetic radiation at the interface between two media with different refractive indices. For ocean surfaces in the infrared, the complex refractive index of seawater is:

```
n̂ = n + ik
```

where `n` is the real part (refraction) and `k` is the imaginary part (absorption).

### 1.2 Implementation in `ocean_emissivity.F90` (lines 459-574)

#### Horizontal (s-polarization, TE) Reflectance:
```fortran
! Fresnel equations for s-wave (TE mode)
num_h = coschi_complex - sqrt_term
den_h = coschi_complex + sqrt_term
rh = abs(num_h/den_h)**2
```

**Mathematical form:**
```
rₕ = |cos(χ) - √(n̂² - sin²χ)| / |cos(χ) + √(n̂² - sin²χ)|
Rₕ = |rₕ|²
```

#### Vertical (p-polarization, TM) Reflectance:
```fortran
! Fresnel equations for p-wave (TM mode)
num_v = refm2*coschi_complex - sqrt_term
den_v = refm2*coschi_complex + sqrt_term
rv = abs(num_v/den_v)**2
```

**Mathematical form:**
```
rᵥ = |n̂² cos(χ) - √(n̂² - sin²χ)| / |n̂² cos(χ) + √(n̂² - sin²χ)|
Rᵥ = |rᵥ|²
```

#### Unpolarized Reflectance:
```fortran
ws%reflections%refle_fresnel(i,j) = 0.5_jprd*(rv + rh)
```

**Mathematical form:**
```
R = (Rᵥ + Rₕ) / 2
```

### 1.3 Key Implementation Details

| Parameter | Value | Location |
|-----------|-------|----------|
| Complex sqrt handling | `sqrt(refm2 - cmplx(sinchi2, 0.0_jprd))` | line 524 |
| Angular grid | 91 zenith × 181 azimuth | N_THETA=91, N_PHI=181 |
| Numerical stability | `coschi_val = max(0.0_jprd, min(1.0_jprd, coschi_val))` | line 509 |
| NaN/Inf check | `ieee_is_finite()` validation | line 550 |

### 1.4 Fresnel Lookup Table (PHASE 2 Optimization)

For performance, a 3D LUT pre-computes Fresnel reflectance:

```fortran
global_fresnel_lut%lut_data(FRESNEL_LUT_SIZE_N, FRESNEL_LUT_SIZE_N, FRESNEL_LUT_SIZE_CHI)
! Dimensions: 256 × 256 × 256 = 16.7 million values
! Range: n_real ∈ [1.0, 3.0], n_imag ∈ [0.0, 1.5], χ ∈ [0, π/2]
```

**Trilinear interpolation** is used for fast lookup (lines 1744-1822).

### 1.5 Polarization Frame Rotation for 2D Surfaces

When integrating over a 2D rough ocean surface, each wave facet can be tilted in the cross-wind direction (φ_n ≠ 0), which rotates the local polarization plane relative to the global (observer-frame) V/H reference. This mixes local V and H contributions when projecting to global coordinates.

**Reference:** Li, Pinel & Bourlier (2012), "Polarized infrared emissivity of 2D sea surfaces with one surface reflection," *Remote Sensing of Environment*, 124, 299-309.

#### Rotation angle β

The rotation angle β between the local and global H-polarization directions is:
```
sin²β = sin²(θ_n) · sin²(φ_n) / sin²(χ)
cos²β = 1 - sin²β
```

where χ is the local incidence angle and (θ_n, φ_n) are the facet normal angles.

#### Rotated emissivities (global frame)

```
ε_V_global = ε_V_local · cos²β + ε_H_local · sin²β
ε_H_global = ε_V_local · sin²β + ε_H_local · cos²β
```

#### Properties
- **Unpolarized invariant:** ε_V_global + ε_H_global = ε_V_local + ε_H_local (cos²β + sin²β = 1)
- **1D surface (φ_n = 0):** sin β = 0, no rotation — reduces to the standard Fresnel result
- **Nadir (θ = 0):** R_V = R_H at normal incidence, so the rotation is moot
- **Effect:** Reduces the V-H polarization contrast, with the largest impact at oblique viewing angles and high wind speeds

#### Implementation in `ocean_emissivity.F90` (lines 545-568)

The rotation is applied in `compute_fresnel_kernel` for first-surface emission (not reflection-angle mode). When sin χ is sufficiently large, sin²β is computed and used to rotate the local Fresnel V/H to global V/H. When sin χ ≈ 0 (near-normal incidence) or in reflection-angle mode, the rotation is skipped since R_V ≈ R_H.

---

## 2. Cox-Munk Integration

### 2.1 Physical Foundation

The Cox-Munk model (1954) describes the statistical distribution of wave facet slopes on a wind-roughened ocean surface. It assumes an **isotropic Gaussian slope distribution**.

### 2.2 Slope Variance (Masuda Constants)

From `ocean_emissivity.F90` (lines 105-106):

```fortran
real(kind=jprd), parameter :: MASUDA_CONST_A = 0.003_jprd
real(kind=jprd), parameter :: MASUDA_CONST_B = 0.00512_jprd
```

**Slope variance formula:**
```
σ²(U₁₀) = 0.003 + 0.00512 × U₁₀
```

where `U₁₀` is wind speed at 10m height (m/s).

### 2.3 Probability Density Function

From `initialize_grid_arrays()` (lines 406-414):

```fortran
if (ws%angles%theta_n(i) >= theta_n_start_deg) then
    ws%angles%PDF(i) = 0.0_jprd
else
    if (ws%angles%costheta_n(i) > 0.0_jprd) then
        ws%angles%PDF(i) = exp(-ws%angles%tan2theta_n(i)/ws%zsig2) / (pi * ws%zsig2)
    else
        ws%angles%PDF(i) = 0.0_jprd
    end if
end if
```

**Mathematical form (isotropic Gaussian):**
```
P(θₙ) = exp(-tan²(θₙ) / σ²) / (π × σ²)
```

### 2.4 Integration Kernel

From `compute_integral_kernel()` (lines 580-622):

```fortran
ws%reflections%integral_fun(i,j) = pdf_val * ws%reflections%coschi(i,j) * costheta4_inv
```

**Mathematical form:**
```
dI(θₙ, φₙ) = P(θₙ) × cos(χ) / cos⁴(θₙ) × dθₙ × dφₙ
```

where:
- `χ` = local incidence angle between viewing ray and facet normal
- `cos(χ) = cos(θ) cos(θₙ) + sin(θ) sin(θₙ) cos(φₙ)` (line 505-506)

### 2.5 Double Integration

From `compute_emissivity_integrals()` (lines 921-993):

```fortran
! Step 1: Integrate over azimuth (φ) for each zenith (θ)
call get_trapz(deg2rad*ws%angles%phi_n, &
               ws%reflections%emiss_fresnel(i,:)*ws%reflections%integral_fun(i,:), &
               inner_sum(i), trapz_error)

! Step 2: Integrate over zenith (θ)
call get_trapz(ws%angles%costheta_n, inner_sum, mean_emiss, trapz_error)

! Normalization factor
mean_emiss = (2.0_jprd/cosvza)*mean_emiss
```

**Mathematical form:**
```
ε = (2/cos θ) ∫∫ ε_fresnel(θₙ,φₙ) × dI(θₙ,φₙ) dθₙ dφₙ
```

Integration limits: θₙ ∈ [0°, 90°], φₙ ∈ [0°, 180°]

---

## 3. Multiple Reflection (Masuda 2006)

### 3.1 Physical Foundation

The Masuda (2006) multiple reflection model accounts for radiation that bounces between adjacent wave facets before escaping. This is significant for:
- High wind speeds (more tilted facets)
- Large viewing angles (more inter-facet interactions)
- High absorption (low wavenumbers where k is large)

### 3.2 Implementation Parameters

From `ocean_emissivity.F90` (lines 81, 117-118):

```fortran
integer(kind=jpim), parameter :: MAX_REFLECTION_ORDER = 5
real(kind=jprd), parameter :: REFLECTION_CONVERGENCE_TOL = 1.0e-8_jprd
integer(kind=jpim), parameter :: MIN_REFLECTION_ORDER = 2
```

### 3.3 Algorithm (`compute_multiple_reflection_kernel`, lines 1085-1203)

#### Iterative Loop:
```fortran
do k = 1, MAX_REFLECTION_ORDER
    ! Compute contribution from k-th reflection order
    ws%iterations%temp_integral(i,j) = ws%iterations%weight(i,j) * &
                                       refle_fresnel0(i,j) * &
                                       ws%iterations%normalized_emiss_sp(i,j) * &
                                       integral_fun0(i,j)

    ! Double integration (φ then θ)
    call get_trapz(deg2rad*ws%angles%phi_n, ws%iterations%temp_integral(i,:), inner_sum(i))
    call get_trapz(ws%angles%costheta_n, inner_sum, refle_od)

    ! Accumulate contribution
    reflection_total = reflection_total + refle_od / norm_factor0

    ! EARLY TERMINATION (PHASE 1 optimization)
    if (abs(refle_od_norm) < REFLECTION_CONVERGENCE_TOL) exit
end do
```

#### Mathematical Series:
```
ε_total = ε₁ + Σₙ₌₂⁵ Rⁿ⁻¹ × ε₁

Where:
ε₁ = first-order (no reflection) emissivity
R = Fresnel reflectance at facet
n = reflection order (up to 5)
```

#### Convergence Check:
```fortran
if (USE_EARLY_TERMINATION .and. k >= MIN_REFLECTION_ORDER) then
    if (ws%last_reflection_error < REFLECTION_CONVERGENCE_TOL) exit
end if
```

### 3.4 Reflection Angle Geometry

From `compute_reflection_angles_kernel()` (lines 776-831):

```fortran
! Reflection direction after bouncing off facet
cosval = -2.0_jprd * ws%angles%costheta_n(i) * ws%reflections%coschi(i,j) + cosvza
ws%reflections%theta_r(i,j) = acos(cosval) * (180.0_jprd/pi)
```

**Physical interpretation:**
The reflected ray direction is computed using the law of reflection:
```
θᵣ = arccos(-2 cos(θₙ) cos(χ) + cos(θ))
```

---

## 4. Seawater Optical Constants

### 4.1 Physical Foundation

The complex refractive index of seawater depends on:
- **Wavenumber** (ν): 10-5000 cm⁻¹
- **Temperature** (T): 273-313 K
- **Salinity** (S): 0-45 PSU

### 4.2 Data Source

From `sea_water_nk.F90` (lines 24-29):

```fortran
! Data Source:
!   - water_optical_constants.nc: 6496 wavenumbers, 11 temperatures
!   - Temperature range: 273-313 K
!   - Salinity range: 0-45 PSU
```

### 4.3 Salinity Correction (Pinkley & Williams, 1976)

From `compute_nk_internal()` (lines 319-326):

```fortran
! Step 1: Apply salinity correction to pure water values
do it = 1, data%ntemp
    refra_real_base(:, it) = data%refra_real_part(it, :) + data%n_salt_corr(:) * salinity_in
    refra_imag_base(:, it) = data%refra_imag_part(it, :) + data%k_salt_corr(:) * salinity_in
end do
```

**Mathematical form:**
```
n(ν,T,S) = n_pure(ν,T) + (dn/dS)(ν) × S
k(ν,T,S) = k_pure(ν,T) + (dk/dS)(ν) × S
```

### 4.4 Log-Space Interpolation

From lines 329-343:

```fortran
! Pre-compute log of wavenumber arrays (log-space interpolation for stability)
wavenumber_log = log(data%wavenumber)
wavenumber_target_log = log(wavenumber_target)

! Step 2: Interpolate over wavenumber for each temperature
do it = 1, data%ntemp
    ! Interpolate real part in log-space
    call get_interpolation_1d(wavenumber_log, log(refra_real_base(:, it)), &
        wavenumber_target_log, temp_real_out_log(:, it), method_interp)
    temp_real_out(:, it) = exp(temp_real_out_log(:, it))

    ! Interpolate imaginary part in log-space
    call get_interpolation_1d(wavenumber_log, log(refra_imag_base(:, it)), &
        wavenumber_target_log, temp_imag_out_log(:, it), method_interp)
    temp_imag_out(:, it) = exp(temp_imag_out_log(:, it))
end do
```

**Why log-space?**
The imaginary part `k` spans 4+ orders of magnitude (0.01-1.0), making linear interpolation numerically unstable. Log-space interpolation preserves accuracy across this wide range.

### 4.5 Temperature Selection

From lines 347-349:

```fortran
! Step 3: Select closest temperature
pos_close_id = minloc(abs(data%temperature - temperature_in), 1)
real_out = temp_real_out(:, pos_close_id)
imag_out = temp_imag_out(:, pos_close_id)
```

**Closest-temperature approach** is used rather than interpolation between temperature grid points.

### 4.6 Spectral Regions

| Region | Wavenumber (cm⁻¹) | n range | k range | Notes |
|--------|------------------|---------|---------|-------|
| Far-IR | 10-200 | 1.5-2.5 | 0.4-1.1 | High absorption |
| Thermal window | 800-1200 | 1.2-1.3 | 0.01-0.1 | Low absorption |
| Near-IR | 2000-5000 | 1.3-1.5 | 0.01-0.4 | Moderate absorption |

---

## 5. Governing Equations Summary

### Ocean Emissivity:
```
ε(θ,ν) = 1 - ∫∫ R(θᵢ,ν) · P(θₙ,φₙ,σ²) · cos(χ) dθₙ dφₙ
```

### Fresnel Reflectance (complex n + ik):
```
rₕ = (n₀ cos θᵢ - n cos θₜ) / (n₀ cos θᵢ + n cos θₜ)  [horizontal, s-pol, TE]
rᵥ = (n cos θᵢ - n₀ cos θₜ) / (n cos θᵢ + n₀ cos θₜ)  [vertical, p-pol, TM]
R = (|rᵥ|² + |rₕ|²) / 2  [unpolarized]
```

### Seawater Refractive Index:
```
n(ν,T,S) = n_pure(ν,T) + (dn/dS)×S
k(ν,T,S) = k_pure(ν,T) + (dk/dS)×S
```

### Planck Function (narrow-band integration):
```
B(ν,T) = c₁ν³ / (exp(c₂ν/T) - 1)
ε̄ = ∫ε(ν)·B(ν,T)dν / ∫B(ν,T)dν
```

---

## 6. Effective View Angle Processing

### 6.1 Physical Motivation

Ocean surface emissivity depends on the viewing geometry relative to wave facets, not on the satellite viewing angle alone. The **effective view angle** (Nalli et al. 2008, 2023) is the single equivalent incidence angle that, when used in a flat-surface Fresnel calculation, reproduces the emissivity of the full rough-surface integration at a given viewing angle and wind speed.

This approach decouples the expensive Cox-Munk integration (Section 2) from operational use: the integration is performed once at high resolution to build a lookup table (LUT), and subsequent calculations use the LUT to convert any (viewing angle, wind speed) pair to an effective angle for a fast Fresnel evaluation.

### 6.2 LUT Generation (Offline)

The effective view angle LUT is generated by the Python script `docs/external_info/process_effective_view_angle.py` using the **Spectral Variance Minimization (SVM)** technique from Nalli et al. (2008) Equation 30 / Nalli et al. (2023) Equation 27:

```
sigma_TB = sqrt( (1/(n-1)) * sum_i [TB_calc(v_i) - TB_ref(v_i)]^2 )
```

**Algorithm:**

1. **Reference data**: Full-physics emissivity `epsilon_ref(v, theta, U)` computed at 76 viewing angles (0-75 deg), 19 wind speeds (0-18 m/s), and 21 temperatures (271-311 K) with complete Cox-Munk integration and multiple reflections.

2. **Calculated data**: Flat-surface Fresnel emissivity `epsilon_calc(v, theta)` computed at 851 angles (0-85 deg in 0.1 deg steps) without wind effects.

3. **Spectral matching**: For each (reference angle, wind speed) pair, search the calculated angle grid to find the angle `theta_eff` that minimizes `sigma_TB` over the atmospheric window [980, 1250] cm^-1. This uses a two-phase search (coarse then fine) with 0.1 deg precision and 0.4% relative bias tolerance.

4. **Multi-temperature processing**: LUTs are generated per temperature, then combined by the merging algorithm (`docs/external_info/update_effective_angle.md`) which consolidates per-temperature results into the full angular range (0-75 deg) using cubic spline transitions and physical constraint enforcement.

The resulting `effective_view_angle.nc` contains a 3D LUT of shape (19 winds, 71 angles, 21 SSTs) mapping geometric viewing angles 0-70 deg to effective viewing angles for each wind speed and sea surface temperature. At runtime, the 3D LUT is interpolated along SST (PCHIP) to produce a 2D (angle, wind) slice. This LUT is wavenumber-independent: the SVM matching in the atmospheric window produces effective angles that generalize across the full infrared spectrum.

### 6.3 LUT Characteristics

The high-resolution 3D LUT (71 angles x 19 winds x 21 SSTs) resolves the effective view angle across the full parameter space with 1-degree angular resolution, 1 m/s wind resolution, and 2 K SST resolution. At runtime, PCHIP interpolation along SST (`interpolate_lut_sst`) produces a 2D (angle, wind) slice tailored to the specific sea surface temperature, capturing the SST-dependent variation of effective view angle that was previously averaged out.

### 6.4 Two-Phase PCHIP LUT Smoothing (`smooth_effective_angle_lut`)

**Location:** `lib/mathlib.F90`

The high-resolution 3D LUT (71 angles x 19 winds x 21 SSTs) resolves effective view angles across the full angular, wind speed, and SST parameter space. After SST interpolation by `interpolate_lut_sst`, the resulting 2D slice is smoothed at runtime using two-phase PCHIP processing applied independently to each wind speed column:

#### Phase 1: PCHIP Re-interpolation

Monotone cubic Hermite (PCHIP) spline interpolation through the 5-degree reference anchors (0, 5, 10, ..., 70 deg) onto the 1-degree grid. This removes sub-anchor sampling noise while preserving the exact LUT values at all anchor points. Fritsch-Carlson guarantees monotonicity between anchors.

```
x_anchor = [0, 5, 10, ..., 65, 70]     (15 points)
y_anchor = eff_angle_LUT(x_anchor, wind) (LUT values)
y_smooth = PCHIP(x_anchor, y_anchor, [0,1,...,70])
```

#### Phase 2: PCHIP Re-fit with Monotonicity Enforcement

Phase 2 guarantees smooth monotonicity by re-fitting through the same 5-degree anchor grid:

1. **Sample** the PCHIP output at the 5-degree anchor points (0, 5, 10, ..., 70)
2. **Enforce** non-negativity and monotonicity on these ~15 sparse anchor values
3. **PCHIP re-interpolate** from the monotone anchors back to the 1-degree grid

The Fritsch-Carlson PCHIP on monotone input guarantees C1-continuous,
monotonically non-decreasing output. Forward-clamping on the sparse anchors
(rather than the dense 1-degree grid) avoids zero-derivative plateaus.

Because the high-resolution LUT fully resolves the effective angle behavior up to 70 deg (including the near-70-degree region), no additional correction is needed at the 70-degree boundary.

### 6.5 Angle Blending at the 70-Degree Boundary

**Location:** `src/driver.F90`

For user-requested angles spanning the 70-degree boundary, the driver splits the computation:

1. **Effective-angle group** (theta <= 70 deg): Emissivity computed using the smoothed effective view angle LUT.

2. **Normal-angle group** (theta >= 70 deg): Emissivity computed using the normal (geometric) view angle LUT from `normal_view_angle.nc`.

3. **Continuity correction** (`blend_angle_results`): A constant offset `diff_70` is computed per wind speed and wavenumber:
   ```
   diff_70 = emiss_eff(70 deg) - emiss_nor(70 deg)
   emiss_final(>70 deg) = emiss_nor(>70 deg) + diff_70
   ```
   With the high-resolution LUT, `diff_70` is small, and the offset correction ensures smooth continuity across the boundary.

4. **Monotonicity enforcement** (`enforce_emissivity_monotonicity`): A final forward-clamping pass ensures emissivity decreases monotonically with increasing angle across the full range. This safety net catches any residual non-monotonicity from numerical noise, primarily at low angles (0-30 deg) where emissivity changes are very small.

   **Polarized mode note:** Monotonicity enforcement is applied only to **unpolarized** and **H-pol** (s-pol, TE) emissivity. It is **not** applied to V-pol (p-pol, TM), because V-pol emissivity physically increases from nadir toward Brewster's angle (~53° for water) before decreasing at large angles. Enforcing monotonic decrease on V-pol would clamp the physical Brewster maximum, producing an incorrect flat curve.

### 6.6 End-to-End Data Flow

```
[Offline: Python]                    [Runtime: Fortran]

SVM spectral matching                effective_view_angle.nc (3D)
(Nalli et al. 2008/2023)   --->     loaded by driver.F90
                                           |
                                    interpolate_lut_sst()
                                     PCHIP interpolation along SST
                                     3D (71 x 19 x 21) --> 2D (71 x 19)
                                           |
                                    smooth_effective_angle_lut()
                                     Phase 1: PCHIP through 5-deg anchors
                                     Phase 2: Non-negativity, monotonicity
                                           |
                                    save_effective_view_angle()
                                     --> effective_view_angle.nc (smoothed LUT)
                                           |
                                    effective_view_angle()
                                     Trilinear interpolation for
                                     user-requested (angle, wind) pairs
                                           |
                              +------------+------------+
                              |                         |
                      theta <= 70 deg            theta > 70 deg
                      Fresnel(eff_angle)         Fresnel(geometric)
                              |                         |
                              +----> blend_angle_results()
                                     diff_70 continuity correction
                                           |
                                    enforce_emissivity_monotonicity()
                                     Forward-clamp safety net
                                     (unpolarized and H-pol only;
                                      V-pol skipped — Brewster maximum)
                                           |
                                    Final emissivity output
```

---

## References

1. **Wu, X., & Smith, W. L. (1997)**. *Emissivity of rough sea surface for 8–13 µm*. Applied Optics, 36(12), 2609-2619.

2. **Masuda, K. (2006)**. *Infrared sea surface emissivity including multiple reflection effect*. Remote Sensing of Environment, 103(4), 488-496.

3. **Cox, C., & Munk, W. (1954)**. *Measurement of the roughness of the sea surface from photographs of the sun's glitter*. Journal of the Optical Society of America, 44(11), 838-850.

4. **Pinkley, L. W., & Williams, D. (1976)**. *Optical properties of sea water in the infrared*. Journal of the Optical Society of America, 66(6), 554-558.

5. **Henderson, B. G., et al. (2003)**. *The polarized emissivity of a wind-roughened sea surface*. Remote Sensing of Environment, 88(4), 453-467.

6. **Nalli, N. R., et al. (2008)**. *Emissivity and reflection model for calculating unpolarized isotropic water surface-leaving radiance*. Applied Optics, 47(21), 3701-3721.

7. **Nalli, N. R., et al. (2023)**. *Reducing biases in thermal infrared surface radiance calculations over global oceans*. IEEE TGRS, 61, 1-18.

8. **Li, Z., Pinel, N., & Bourlier, C. (2012)**. *Polarized infrared emissivity of 2D sea surfaces with one surface reflection*. Remote Sensing of Environment, 124, 299-309.

9. **Wang, S., Yang, P., Brindley, H. E., Huang, X., & L'Ecuyer, T. S. (2026)**. *Enhanced full spectral temperature-dependent refractive index of liquid water from supercooled to ambient conditions*. Geophysical Research Letters, 53(6), e2025GL119385. https://doi.org/10.1029/2025GL119385
