# IOSEE Architecture Documentation

This directory contains documentation of the IOSEE codebase architecture, performance optimizations, and implementation patterns.

## Contents

- [**architecture_overview.md**](architecture_overview.md) - Comprehensive codebase analysis:
  - Directory structure and module dependencies
  - Key data structures
  - Performance optimization strategies
  - Build system and compiler support

## Quick Statistics

| Metric | Value |
|--------|-------|
| Fortran Lines | 11,018 |
| Data Files | 127 MB |
| Performance | 12,000+ calc/sec |

## Module Structure

```
utils.F90 → configuration.F90 → mathlib.F90
                                    ↓
netcdf_handler.F90 ─┬─→ sea_water_nk.F90
                    └─→ ocean_emissivity.F90
                              ↓
                        driver.F90 → save_output.F90
```

## Performance Optimizations

| Phase | Technique | Speedup |
|-------|-----------|---------|
| 1 | Adaptive grid | 4X |
| 2 | Fresnel LUT | 2X |
| 3 | Early termination | 1.5X |

## See Also

- [Physics Deep-Dive](../physics/physics_deep_dive.md)
- [User Guide](../user_guide/)
- [API Documentation](../api/)
