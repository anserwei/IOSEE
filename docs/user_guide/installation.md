# IOSEE Installation Guide

IOSEE is a standalone Fortran executable built with CMake.

---

# Installation

## Prerequisites

| Software | Minimum Version | Purpose |
|----------|-----------------|---------|
| CMake | 3.16+ | Build system |
| Fortran Compiler | gfortran 12+ / ifort 19+ / ifx 2024+ | Compilation |
| NetCDF-Fortran | 4.5+ | Data I/O |
| OpenMP | 4.0+ | Parallelization (optional) |

## Install Dependencies

### macOS (Homebrew)
```bash
brew install gcc netcdf cmake
```

### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install gfortran libnetcdf-dev libnetcdff-dev cmake
```

### CentOS/Rocky Linux
```bash
sudo yum install gcc-gfortran netcdf-fortran-devel cmake
```

### HPC Systems (module environment)
```bash
module load gcc netcdf cmake
# Or: module load intel netcdf cmake
```

## Build

```bash
# Clone or download IOSEE
cd /path/to/IOSEE

# Build release version (recommended)
make release

# Verify build
./examples/iosee examples/spectral_option1.config
```

### Build Options

| Target | Description |
|--------|-------------|
| `make release` | Optimized build with OpenMP (recommended) |
| `make debug` | Debug symbols, bounds checking |
| `make serial` | No OpenMP (single-threaded) |
| `make reproducible` | Bit-for-bit reproducible results |
| `make clean` | Remove build artifacts |

### Compiler Selection

```bash
# GNU Fortran (default)
FC=gfortran make release

# Intel Classic Fortran
FC=ifort make release

# Intel LLVM Fortran
FC=ifx make release

# Cray Fortran (HPC)
FC=ftn make release
```

## Run

```bash
cd examples

# Run with configuration file
./iosee spectral_option1.config

# Specify output file
./iosee spectral_option1.config output=my_output.nc
```

Running `./iosee` with no arguments prints the usage banner.

### Example Output
```
IOSEE v1.0.0 - Ocean Surface Emissivity Calculator
Reading configuration: spectral_option1.config
Computing spectral emissivity...
  Wavenumber range: 800.0 - 1200.0 cm^-1 (401 points)
  Temperature: 300.0 K, Salinity: 35.0 PSU
  Angles: 4, Wind speeds: 4
Output written to: ocean_emissivity_spectral.nc
```

## Verify Installation

```bash
cd examples

# Run a test case
./iosee spectral_option1.config

# Check output was created
ls -la ocean_emissivity_spectral.nc

# Optionally inspect with ncdump
ncdump -h ocean_emissivity_spectral.nc
```

---

# Troubleshooting

## Fortran Mode Issues

### "NetCDF not found"

```bash
# Check if nf-config exists
which nf-config

# Set path manually if needed
export NETCDF_FORTRAN_DIR=/path/to/netcdf-fortran
make release
```

### "OpenMP not supported"

```bash
# Build without OpenMP
make serial
```

### Apple Silicon (M1/M2/M3/M4)

```bash
# Use Homebrew gfortran
brew install gcc
FC=/opt/homebrew/bin/gfortran make release
```

---

# Uninstallation

```bash
cd /path/to/IOSEE
make clean
```

---

**Next:** See [Configuration Guide](configuration.md) for configuration file options.
