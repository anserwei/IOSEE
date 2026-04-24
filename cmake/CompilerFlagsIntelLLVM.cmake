#==============================================================================
# Intel LLVM Fortran Compiler Configuration (ifx)
# For Intel ifx (2024.0+) with LLVM backend
# CMake identifies this as CMAKE_Fortran_COMPILER_ID = "IntelLLVM"
#==============================================================================

message(STATUS "-- Configuring Intel ifx (LLVM) compiler flags")

# Basic compiler flags with platform-specific preprocessor definitions
set(FORTRAN_FLAGS_BASE "-fpp -assume realloc_lhs -heap-arrays 64")

# Add platform-specific preprocessor flags for proper platform detection
if(UNIX AND NOT APPLE)
    # Linux/Unix systems (including HPC)
    set(FORTRAN_FLAGS_BASE "${FORTRAN_FLAGS_BASE} -D__linux__ -Dlinux")
elseif(APPLE)
    # Apple systems
    set(FORTRAN_FLAGS_BASE "${FORTRAN_FLAGS_BASE} -D__APPLE__")
endif()

# Debug configuration
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -check all -traceback -warn all ${FORTRAN_FLAGS_BASE}")

# Release flags - platform-aware optimization
if(APPLE)
    # macOS: use host-specific optimizations
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost ${FORTRAN_FLAGS_BASE}")
elseif(HPC_ENVIRONMENT)
    # HPC environments: use conservative AVX2 for broad compatibility
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -xCORE-AVX2 ${FORTRAN_FLAGS_BASE}")
else()
    # Standard Linux: use native optimizations
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost ${FORTRAN_FLAGS_BASE}")
endif()

# Reproducibility mode - strict IEEE 754 compliance
if(ENABLE_REPRODUCIBLE)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fp-model=strict -fimf-arch-consistency=true")
    message(STATUS "-- Reproducibility: ENABLED (strict IEEE 754)")
else()
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -fp-model=fast")
endif()

# OpenMP parallelization
if(ENABLE_OPENMP)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -qopenmp")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -qopenmp")
endif()

# Link-Time Optimization (LTO)
if(ENABLE_LTO)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -ipo")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ipo")
endif()

# Add NetCDF include path after flags are set
if(NETCDF_INCLUDE)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${NETCDF_INCLUDE}")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${NETCDF_INCLUDE}")
endif()

# Set unified flags for build type
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEBUG}")
else()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_RELEASE}")
endif()
