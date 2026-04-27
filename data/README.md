# IOSEE — `data/` directory

This directory holds the physical lookup tables used by the IOSEE driver at
runtime.

## Files tracked in git

| File                                | Size   | Purpose                                                    |
|-------------------------------------|--------|------------------------------------------------------------|
| `effective_view_angle.nc`           | 84 KB  | Derived effective view angles.                             |
| `normal_view_angle.nc`              | 10 KB  | Geometric view angles.                                     |
| `water_optical_constants_HQ.nc`     | 191 KB | Reduced-precision pure water refractive index (HQ).        |
| `water_optical_constants_HQ_full.nc`| 191 KB | Full-resolution high-quality refractive index (HQ).        |
| `water_optical_constants_HQ_SE.nc`  | 1.3 MB | Pure water refractive index (HQ-SE).                       |

## Large file — download separately

The primary input, `water_optical_constants.nc` (≈ 128 MB), is **not tracked
in git** because it exceeds GitHub's 100 MB per-file limit. It is distributed
as a GitHub Release asset.

Download it into this directory before running IOSEE:

```bash
cd data
curl -LO https://github.com/anserwei/IOSEE/releases/download/v1.0.0/water_optical_constants.nc
# or:
wget    https://github.com/anserwei/IOSEE/releases/download/v1.0.0/water_optical_constants.nc
```

Optionally verify the download:

```bash
shasum -a 256 water_optical_constants.nc
```

## Provenance

`water_optical_constants.nc` derives from:

> Wang, S., Yang, P., Brindley, H. E., Huang, X., & L'Ecuyer, T. S. (2026).
> *Enhanced full spectral temperature-dependent refractive index of liquid
> water from supercooled to ambient conditions.* Geophysical Research Letters,
> 53(6), e2025GL119385. https://doi.org/10.1029/2025GL119385

`effective_view_angle.nc` derives from Nalli et al. (2008, 2023) via
spectral variance minimization and is further smoothed at runtime; see
`docs/physics/physics_deep_dive.md` §6.
