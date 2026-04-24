#==============================================================================
# GNU Fortran Compiler Configuration
#==============================================================================


# Basic compiler flags with platform-specific preprocessor definitions
set(BASIC_FLAGS "-fno-range-check -cpp -fallow-argument-mismatch -std=f2008")

# Add platform-specific preprocessor flags for proper platform detection
if(UNIX AND NOT APPLE)
    # Linux/Unix systems (including HPC)
    set(BASIC_FLAGS "${BASIC_FLAGS} -D__linux__ -Dlinux")
elseif(APPLE)
    # Apple systems
    set(BASIC_FLAGS "${BASIC_FLAGS} -D__APPLE__")
endif()

# Universal Architecture-Specific Optimization
if(APPLE_SILICON)
    # Apple Silicon: use -mcpu instead of -march for ARM64
    set(OPT_FLAGS "-O3 -mcpu=native")
elseif(APPLE_INTEL)
    # Intel Mac: use native optimizations for maximum performance
    set(OPT_FLAGS "-O3 -march=native -mtune=native")
elseif(HPC_ENVIRONMENT)
    # HPC environments: use compatible optimizations across node types
    if(ARCH MATCHES "x86_64")
        # Use conservative x86-64 baseline for HPC cluster compatibility
        set(OPT_FLAGS "-O3 -march=x86-64 -mtune=generic")
    elseif(ARCH MATCHES "aarch64")
        # ARM64 HPC systems (e.g., Fugaku, ARM-based clusters)
        set(OPT_FLAGS "-O3 -march=armv8-a -mtune=generic")
    else()
        # Conservative optimization for other HPC architectures
        set(OPT_FLAGS "-O3")
    endif()
elseif(ARCH MATCHES "x86_64")
    # Standard Linux x86_64: use native optimizations
    set(OPT_FLAGS "-O3 -march=native -mtune=native")
elseif(ARCH MATCHES "aarch64")
    # Standard Linux ARM64: use native optimizations
    set(OPT_FLAGS "-O3 -march=native -mtune=native")
else()
    # Fallback for unknown architectures
    set(OPT_FLAGS "-O3")
endif()

# Advanced vectorization and loop optimizations
set(OPT_FLAGS "${OPT_FLAGS} -ftree-vectorize")
if(NOT HPC_ENVIRONMENT)
    # Aggressive optimizations - may cause issues on some HPC systems
    set(OPT_FLAGS "${OPT_FLAGS} -fvect-cost-model=unlimited")
    set(OPT_FLAGS "${OPT_FLAGS} -funroll-loops -funroll-all-loops")
    set(OPT_FLAGS "${OPT_FLAGS} -floop-nest-optimize -ftree-loop-distribution")
else()
    # Conservative vectorization for HPC compatibility
    set(OPT_FLAGS "${OPT_FLAGS} -funroll-loops")
endif()
set(OPT_FLAGS "${OPT_FLAGS} -ftree-loop-vectorize -ftree-loop-optimize")

# Math and floating point optimizations for scientific computing
if(ENABLE_REPRODUCIBLE)
    # Reproducibility mode: strict IEEE 754 compliance
    set(OPT_FLAGS "${OPT_FLAGS} -ffp-contract=off -fno-fast-math")
    message(STATUS "-- Reproducibility: ENABLED (strict IEEE 754, no fast-math)")
elseif(HPC_ENVIRONMENT)
    # Conservative math optimizations for HPC environments
    set(OPT_FLAGS "${OPT_FLAGS} -fno-math-errno")
else()
    # Aggressive math optimizations for local development
    set(OPT_FLAGS "${OPT_FLAGS} -ffast-math -fno-math-errno -ffinite-math-only")
    set(OPT_FLAGS "${OPT_FLAGS} -fno-signed-zeros -fno-trapping-math")
endif()

# Cache and memory optimizations
set(OPT_FLAGS "${OPT_FLAGS} -fprefetch-loop-arrays -foptimize-strlen")
set(OPT_FLAGS "${OPT_FLAGS} -falign-functions=32 -falign-loops=32")

# Function inlining and interprocedural optimizations
set(OPT_FLAGS "${OPT_FLAGS} -finline-functions -finline-small-functions")
set(OPT_FLAGS "${OPT_FLAGS} -fipa-pta -fdevirtualize-speculatively")

# SIMD optimizations
set(OPT_FLAGS "${OPT_FLAGS} -fsimd-cost-model=unlimited")

# Warning flags (suppress common false positives)
set(WARN_FLAGS "-Wall -Wno-tabs -Wno-unused-label -Wno-unused-dummy-argument")
set(WARN_FLAGS "${WARN_FLAGS} -Wno-unused-variable -Wno-compare-reals -Wno-maybe-uninitialized")

# Debug configuration
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(DEBUG_FLAGS "-fcheck=bounds -fbacktrace -g -fno-inline")
    set(OPT_FLAGS "-O1")  # Reduce optimization for debugging
    set(LTO_FLAGS "")
else()
    set(DEBUG_FLAGS "")
    # Link Time Optimization (LTO) for maximum performance
    if(HPC_ENVIRONMENT)
        # HPC environments: LTO can cause linking issues, disable by default
        set(LTO_FLAGS "")
    elseif(APPLE_SILICON)
        # Apple Silicon: simplified LTO flags
        set(LTO_FLAGS "-flto")
    else()
        # Intel/AMD: full LTO with linker plugin
        set(LTO_FLAGS "-flto=auto -fuse-linker-plugin")
    endif()
endif()

# Combine all flags  
set(CMAKE_Fortran_FLAGS "${BASIC_FLAGS} ${OPT_FLAGS} ${WARN_FLAGS} ${DEBUG_FLAGS} ${OMP_FLAGS}")

# Add NetCDF include path after flags are set
if(NETCDF_INCLUDE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${NETCDF_INCLUDE}")
endif()

# Linker flags
if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LTO_FLAGS}")
endif()
