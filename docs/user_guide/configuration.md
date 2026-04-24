# IOSEE Configuration Guide

## Configuration File Format

IOSEE uses Fortran namelist format for configuration files. All parameters are specified within a `&configuration` block.

## Basic Configuration

```fortran
&configuration
! Mode selection
mode = 'spectral'           ! 'spectral' or 'narrow-band'
polarization = 'disable'    ! 'enable' or 'disable'

! Spectral range
wavenum_start = 10.0        ! Start wavenumber [cm^-1]
wavenum_end = 3250.0        ! End wavenumber [cm^-1]
n_wavenum = 3241            ! Number of wavenumber points

! Environmental conditions
sea_surface_temperature = 300.0  ! SST [K]
sea_surface_salinity = 35.0      ! Salinity [PSU]

! Viewing geometry
view_angle = 0.0, 30.0, 60.0     ! Viewing angles [degrees]

! Wind conditions
wind = 0.0, 5.0, 10.0, 15.0      ! Wind speeds [m/s]

! Data files
refractive_index_filename = 'data/water_optical_constants.nc'
/
```

## Configuration Parameters

### Mode Selection

| Parameter | Values | Description |
|-----------|--------|-------------|
| `mode` | `'spectral'`, `'narrow-band'` | Calculation mode |
| `polarization` | `'enable'`, `'disable'` | V/H polarization output |

### Spectral Range

| Parameter | Type | Units | Description |
|-----------|------|-------|-------------|
| `wavenum_start` | real | cm^-1 | Start wavenumber (or array for narrow-band) |
| `wavenum_end` | real | cm^-1 | End wavenumber (or array for narrow-band) |
| `n_wavenum` | integer | - | Number of spectral points (spectral mode only) |

**Valid range:** 10 - 3250 cm^-1 (corresponding to 3.1 - 1000 μm)

### Environmental Conditions

| Parameter | Type | Units | Valid Range | Description |
|-----------|------|-------|-------------|-------------|
| `sea_surface_temperature` | real | K | 270 - 310 | Sea surface temperature |
| `sea_surface_salinity` | real | PSU | 0 - 40 | Sea surface salinity |

### Viewing Geometry

| Parameter | Type | Units | Valid Range | Description |
|-----------|------|-------|-------------|-------------|
| `view_angle` | real array | degrees | 0 - 85 | Viewing angles from nadir |

**Note:** Up to 200 angles can be specified. Angles <= 70 deg use effective view angles (Nalli et al. 2008, 2023); angles > 70 deg use normal (geometric) view angles with a continuity correction at the 70 deg boundary.

### Wind Conditions

| Parameter | Type | Units | Valid Range | Description |
|-----------|------|-------|-------------|-------------|
| `wind` | real array | m/s | 0 - 18 | Wind speeds at 10m height |

**Note:** Up to 200 wind speeds can be specified.

### Output Options

| Parameter | Values | Description |
|-----------|--------|-------------|
| `no_multiple_reflection` | `.true.`, `.false.` | Skip multiple reflection calculation |
| `use_effective_angle` | `.true.`, `.false.` | Use effective view angle LUT (default: `.true.`) |

**`use_effective_angle`:** When `.true.` (default), angles <= 70 deg use the Nalli et al. effective view angle LUT and angles > 70 deg use the normal (geometric) view angle LUT with continuity correction at the 70 deg boundary. When `.false.`, the effective view angle LUT is skipped entirely and only the normal (geometric) view angle LUT is used for all angles.

## Example Configurations

### Spectral Mode (Full Spectrum)

```fortran
&configuration
mode = 'spectral'
wavenum_start = 10.0
wavenum_end = 3250.0
n_wavenum = 3241
sea_surface_temperature = 298.15
sea_surface_salinity = 35.0
view_angle = 0.0
wind = 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0
refractive_index_filename = 'data/water_optical_constants.nc'
/
```

### Polarized Mode

```fortran
&configuration
mode = 'spectral'
polarization = 'enable'
wavenum_start = 800.0
wavenum_end = 1200.0
n_wavenum = 401
sea_surface_temperature = 300.0
sea_surface_salinity = 35.0
view_angle = 0.0, 15.0, 30.0, 45.0, 55.0, 65.0
wind = 0.0, 5.0, 10.0, 15.0
refractive_index_filename = 'data/water_optical_constants.nc'
/
```

### Narrow-Band Mode (RRTMGP bands)

```fortran
&configuration
mode = 'narrow-band'
wavenum_start = 10, 250, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 2250, 2380, 2600
wavenum_end = 250, 500, 630, 700, 820, 980, 1080, 1180, 1390, 1480, 1800, 2080, 2250, 2380, 2600, 3250
sea_surface_temperature = 288.0
sea_surface_salinity = 35.0
view_angle = 0.0, 53.0
wind = 0.0, 5.0, 10.0
refractive_index_filename = 'data/water_optical_constants.nc'
/
```

## Command Line Usage (Fortran mode)

```bash
# Basic usage
./iosee config_file.config

# Specify output file
./iosee config_file.config output=my_output.nc

# Examples
./iosee spectral_option1.config
./iosee spectral_option1.config output=result.nc
```

## Output Files

IOSEE produces NetCDF files with the following structure:

### Smoothed Effective View Angle LUT

Every run also writes `effective_view_angle.nc` to the current working directory. This contains the PCHIP-smoothed effective view angle LUT used by the calculation:

| Variable | Dimensions | Description |
|----------|------------|-------------|
| `reference_angle` | (71) | Reference viewing angles 0-70 deg |
| `reference_wind_speed` | (19) | Reference wind speeds 0-18 m/s |
| `effective_view_angle` | (71, 19) | Smoothed effective viewing angles [deg] |

### Spectral Mode Output

| Variable | Dimensions | Description |
|----------|------------|-------------|
| `emissivity` | (n_angle, n_winds, n_wavenumber) | Ocean surface emissivity |
| `wavenumber` | (n_wavenumber) | Wavenumber coordinate [cm^-1] |
| `angle` | (n_angle) | Viewing angle coordinate [degrees] |
| `wind_speed` | (n_winds) | Wind speed coordinate [m/s] |
| `sea_surface_temperature` | scalar | SST [K] |
| `sea_surface_salinity` | scalar | Salinity [PSU] |

### Polarized Mode Output

Additional variables:
| Variable | Description |
|----------|-------------|
| `emissivity_v` | Vertical polarization emissivity |
| `emissivity_h` | Horizontal polarization emissivity |

---

**Next:** See [Physics Deep-Dive](../physics/physics_deep_dive.md) for algorithm details.
