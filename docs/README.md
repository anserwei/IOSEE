# IOSEE Documentation

Comprehensive documentation for the Infrared Ocean Surface Effective Emissivity package.

---

## Documentation Sections

### [Physics](physics/)
Detailed technical analysis of the physical models implemented in IOSEE:
- **[Physics Deep-Dive](physics/physics_deep_dive.md)** - Fresnel equations, Cox-Munk integration, multiple reflections, optical constants, effective view angle processing

### [Architecture](architecture/)
Codebase structure and implementation patterns:
- **[Architecture Overview](architecture/architecture_overview.md)** - Module dependencies, data structures, performance optimizations

### [User Guide](user_guide/)
Getting started and configuration:
- **[Installation](user_guide/installation.md)** - Build and install instructions
- **[Configuration](user_guide/configuration.md)** - Config file format and options
- **[Configurable Flags Reference](user_guide/configurable_flag.docx)** - Complete table of every namelist flag with valid options and defaults (Word `.docx`)

### [API Reference](api/)
Auto-generated API documentation:
- **[API README](api/README.md)** - Instructions for generating FORD documentation

### [Configuration Examples](config_case/)
Sample configuration files for different use cases:
- Spectral mode (wavenumber and wavelength)
- Narrow-band mode (RRTMGP bands)
- Polarized mode (V and H components)

---

## Quick Links

| Topic | Document |
|-------|----------|
| Building IOSEE | [Installation Guide](user_guide/installation.md) |
| Config file syntax | [Configuration Guide](user_guide/configuration.md) |
| Fresnel equations | [Physics Deep-Dive](physics/physics_deep_dive.md#1-fresnel-equations-implementation) |
| Cox-Munk model | [Physics Deep-Dive](physics/physics_deep_dive.md#2-cox-munk-integration) |
| Multiple reflections | [Physics Deep-Dive](physics/physics_deep_dive.md#3-multiple-reflection-masuda-2006) |
| Optical constants | [Physics Deep-Dive](physics/physics_deep_dive.md#4-seawater-optical-constants) |
| Effective view angles | [Physics Deep-Dive](physics/physics_deep_dive.md#6-effective-view-angle-processing) |
| Performance | [Architecture Overview](architecture/architecture_overview.md#3-performance-optimizations) |

---

## Document Summary

| Document | Purpose |
|----------|---------|
| [physics/physics_deep_dive.md](physics/physics_deep_dive.md) | Complete physics implementation details |
| [architecture/architecture_overview.md](architecture/architecture_overview.md) | Codebase structure and optimizations |
| [user_guide/installation.md](user_guide/installation.md) | Build and install instructions |
| [user_guide/configuration.md](user_guide/configuration.md) | Config file format |
| [user_guide/configurable_flag.docx](user_guide/configurable_flag.docx) | Complete flag reference (`.docx`) |
| [api/README.md](api/README.md) | API documentation generation |

---

## Key References

1. **Wu, X., & Smith, W. L. (1997)**. *Emissivity of rough sea surface for 8–13 µm*. Applied Optics, 36(12), 2609-2619.

2. **Masuda, K. (2006)**. *Infrared sea surface emissivity including multiple reflection effect*. Remote Sensing of Environment, 103(4), 488-496.

3. **Cox, C., & Munk, W. (1954)**. *Measurement of the roughness of the sea surface from photographs of the sun's glitter*. JOSA, 44(11), 838-850.

4. **Nalli, N. R., et al. (2023)**. *Reducing biases in thermal infrared surface radiance calculations over global oceans*. IEEE TGRS, 61, 1-18.

---

**Version**: 1.0.0
**License**: Apache License 2.0
