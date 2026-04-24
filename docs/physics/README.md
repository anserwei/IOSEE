# IOSEE Physics Documentation

This directory contains detailed documentation of the physical models and algorithms implemented in IOSEE.

## Contents

- [**physics_deep_dive.md**](physics_deep_dive.md) - Comprehensive analysis of all physical models:
  - Fresnel equations implementation
  - Cox-Munk wave slope integration
  - Multiple reflection (Masuda 2006) algorithm
  - Seawater optical constants processing
  - Effective view angle processing (LUT generation, smoothing, boundary blending)

## Physical Models

| Model | Reference | Description |
|-------|-----------|-------------|
| Ocean Emissivity | Wu & Smith (1997) | Rough sea surface emissivity |
| Multiple Reflections | Masuda (2006) | Inter-facet reflection corrections |
| Wave Slopes | Cox-Munk (1954) | Gaussian slope distribution |
| Optical Constants | Pinkley & Williams (1976) | Salinity-corrected n,k values |
| Effective View Angle | Nalli et al. (2008, 2023) | SST-dependent 3D LUT with PCHIP smoothing |

## Key Equations

### Fresnel Reflectance
```
R = (|rᵥ|² + |rₕ|²) / 2
```

### Cox-Munk Slope Variance
```
σ² = 0.003 + 0.00512 × U₁₀
```

### Seawater Refractive Index
```
n(ν,T,S) = n_pure(ν,T) + (dn/dS) × S
```

## See Also

- [Architecture Overview](../architecture/architecture_overview.md)
- [User Guide](../user_guide/)
- [API Documentation](../api/)
