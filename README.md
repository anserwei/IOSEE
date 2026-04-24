# Infrared Ocean Surface Effective Emissivity package (IOSEE)

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

A high-performance Fortran package for calculating infrared ocean surface emissivity based on **Masuda (2006)** physical models. Supports spectral, narrow-band, and polarized emissivity calculations with OpenMP parallelization.

---

**Version**: 1.0.0
**License**: Apache License 2.0
**Copyright**: (C) 2025 - Texas A&M University
**Maintainer**: Atmospheric & Oceanic Optics Group, Department of Atmospheric Sciences
**Contact**: Jian Wei ; Email: anser@tamu.edu

## Reference

- **Wei, J., Yang, P. (2026)**. *IOSEE: An Efficient and Flexible Infrared Ocean Surface Effective Emissivity Package for Remote Sensing, Numerical Weather Prediction, and Climate Modeling Applications *. in praperation

## Key Features

- **Three calculation modes**: Spectral, narrow-band, and polarized emissivity
- **Full spectral range**: 10-5000 cm^-1 (far-infrared to thermal infrared)
- **High performance**: 12,000+ calculations/second with OpenMP optimization
- **Physical base**: Surface emission, multiple surface reflection based ocean surface emissivity physical models with effective viewing zenith angles approach
- **Scientific accuracy**: Validation against independent shipborne and ground-based spectroradiometric observations spanning diverse oceanic thermal dynamical conditions
- **Cross-platform**: Tested on Apple Silicon

---

## Importance Notes
- **The default water optical constant file**: the water_optical_constants.nc (≈ 128 MB), is not tracked in git; it is located at https://github.com/anserwei/IOSEE/releases/download/v1.0.0/
- **link**: wget https://github.com/anserwei/IOSEE/releases/download/v1.0.0/water_optical_constants.nc

---

## Quick Start

### Prerequisites

| Software | Version | Installation |
|----------|---------|--------------|
| Fortran Compiler | gfortran 12+ / ifort 19+ / ifx 2024+ | System package manager |
| NetCDF-Fortran | 4.5+ | `brew install netcdf` / `apt install libnetcdff-dev` |
| CMake | 3.16+ | `brew install cmake` / `apt install cmake` |
| OpenMP | 4.0+ | Included with compiler (optional) |

### Build

```bash
# Install dependencies (choose your platform)
# macOS:   brew install gcc netcdf cmake
# Ubuntu:  sudo apt install gfortran libnetcdf-dev libnetcdff-dev cmake
# Conda:   conda install -c conda-forge gfortran netcdf-fortran cmake

# Build IOSEE
make release

# Verify build
./examples/iosee examples/spectral_option1.config
```

### Run

```bash
cd examples

./iosee <user_defined>.config output=<user_defined>.nc

```

### Configuration

Create a `.config` file with Fortran namelist syntax:

```fortran
&configuration
mode = 'spectral'
wavenum_start = 800.0
wavenum_end = 1200.0
n_wavenum = 401
sea_surface_temperature = 300.0
sea_surface_salinity = 35.0
view_angle = 0.0, 30.0, 45.0, 60.0
wind = 0.0, 5.0, 10.0, 15.0
/
```

See [docs/user_guide/configuration.md](docs/user_guide/configuration.md) for all options.

---

## File Structure

```
IOSEE/
├── src/                        # Core Fortran source
│   ├── driver.F90              # Main program
│   ├── ocean_emissivity.F90    # Emissivity physics engine
│   └── sea_water_nk.F90        # Optical constants
├── lib/                        # Support libraries
│   ├── mathlib.F90             # Numerical algorithms
│   ├── configuration.F90       # Config parsing
│   ├── netcdf_handler.F90      # NetCDF I/O
│   └── utils.F90               # Utilities
├── data/                       # Optical constants data
├── examples/                   # Example configs
└── docs/                       # Documentation
```

---

## Documentation

| Document | Description |
|----------|-------------|
| [Installation Guide](docs/user_guide/installation.md) | Detailed installation instructions |
| [Configuration Guide](docs/user_guide/configuration.md) | Config file options |
| [Configurable Flags](docs/user_guide/configurable_flag.docx) | Complete namelist flag reference (`.docx`) |
| [Physics Deep-Dive](docs/physics/physics_deep_dive.md) | Algorithm details |
| [Changelog](CHANGELOG.md) | Version history |

---

## Compiler Support

| Compiler | CMake ID | Versions | Reproducibility | OpenMP |
|----------|----------|----------|-----------------|--------|
| GNU gfortran | GNU | 12.0+ | `-ffp-contract=off` | `-fopenmp` |
| Intel ifort | Intel | 19.0+ | `-fp-model strict` | `-qopenmp` |
| Intel ifx | IntelLLVM | 2024.0+ | `-fp-model=strict` | `-qopenmp` |
| Cray ftn | Cray | CCE 15+ | `-h fp0` | `-h omp` |

## Platform Support

| Platform | Architecture | Optimization |
|----------|--------------|--------------|
| macOS | Apple Silicon M1-M4 | `-mcpu=native` |
| macOS | Intel x86_64 | `-march=native` |
| Linux | x86_64 | `-march=native` |
| Linux | ARM64/aarch64 | `-mcpu=native` |
| HPC | Various | Conservative `-O2` |



