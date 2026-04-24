# IOSEE TODO List

This document tracks planned improvements and future development items for the
IOSEE (Infrared Ocean Surface Effective Emissivity)
package.

---

## High Priority

### Documentation
- [ ] **Generate API documentation with FORD**
  - [ ] Add FORD-compatible docstrings to all public subroutines
  - [ ] Generate HTML documentation
  - [ ] Host on GitHub Pages or Read the Docs

---

## Medium Priority

### Code Quality
- [ ] **Improve error handling system**
  - [ ] Add optional logging system with verbosity levels

### Performance Optimizations
- [ ] **GPU acceleration (OpenACC/CUDA)**
  - [ ] Evaluate feasibility for Fresnel LUT interpolation
  - [ ] Consider NVIDIA CUDA Fortran for HPC systems
  - [ ] Benchmark against CPU-only OpenMP version

### User Experience
- [ ] **Add command-line interface improvements**
  - [ ] Add `--help` and `--version` flags
  - [ ] Add progress bar for long calculations

---

## Low Priority

### Longterm Possible Features (Hard)
- [ ] **Add foam emissivity model**
  - [ ] Implement Koepke (1984) foam coverage model
  - [ ] Add whitecap fraction as function of wind speed
  - [ ] Blend foam and water emissivity

- [ ] **Add sea ice emissivity option**
  - [ ] Implement first-year ice optical properties
  - [ ] Add multi-year ice support
  - [ ] Add ice fraction blending with open water

- [ ] **Add atmospheric path correction**
  - [ ] Implement simple atmospheric transmission model
  - [ ] Add water vapor absorption effects
  - [ ] Interface with LBLRTM or RRTMGP

### Distribution
- [ ] **Create Docker container**
  - [ ] Dockerfile with all dependencies
  - [ ] Pre-built images for CI/CD
  - [ ] Singularity container for HPC


---

## Contributing

To contribute to IOSEE development:

1. Check this TODO list for open items
2. Open an issue to discuss proposed changes
3. Fork the repository and create a feature branch
4. Submit a pull request with tests and documentation

---

**Last Updated:** 2026-03-03
**Maintainer:** Dr. Jian Wei (anser@tamu.edu)
**Institution:** Atmospheric & Oceanic Optics Group, Department of Atmospheric Sciences, Texas A&M University
