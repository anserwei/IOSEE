# IOSEE API Documentation

This directory contains auto-generated API documentation for IOSEE.

## Generating Documentation

To generate the API documentation, install FORD and run:

```bash
# Install FORD (https://github.com/Fortran-FOSS-Programmers/ford)
ford ford.md

# Open in browser
open docs/api/index.html
```

## Documentation Contents

After generation, this directory will contain:

```
docs/api/
├── index.html           # Main documentation page
├── lists/
│   ├── modules.html     # List of all modules
│   ├── procedures.html  # List of all subroutines/functions
│   └── types.html       # List of all derived types
├── module/
│   ├── ocean_emissivity.html
│   ├── sea_water_nk.html
│   ├── mathlib.html
│   ├── configuration.html
│   ├── netcdf_handler.html
│   ├── save_output.html
│   └── utils.html
├── proc/
│   ├── get_emissivity_optimized.html
│   ├── get_polarized_emissivity_optimized.html
│   └── ... (all public subroutines)
├── type/
│   ├── config_type.html
│   ├── thread_workspace.html
│   └── error_type.html
└── search.html          # Search interface
```

## FORD Docstring Format

IOSEE uses FORD-compatible docstrings:

```fortran
!> Calculate ocean surface emissivity using Wu & Smith (1997) model.
!>
!> This subroutine computes spectral emissivity for a wind-roughened
!> ocean surface including multiple reflection effects.
!>
!> @param refm Complex refractive index of sea water
!> @param effective_angles_2d Effective viewing angles [degrees]
!> @param pwind Wind speed at 10m [m/s]
!> @param emissivity_final Output emissivity array [0-1]
!>
!> @note Requires prior initialization via initialize_workspace_pool_emiss
!>
!> @see Masuda, K. (2006). Infrared sea surface emissivity including
!>      multiple reflection effect. RSE, 103(4), 488-496.
subroutine get_emissivity_optimized(refm, effective_angles_2d, pwind, emissivity_final)
```

---

**Note:** This README will be replaced when you run `ford ford.md`.
