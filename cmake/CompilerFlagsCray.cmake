#==============================================================================
# Cray Fortran Compiler Configuration (ftn)
# For Cray Compiling Environment (CCE) 15.0+
# CMake identifies this as CMAKE_Fortran_COMPILER_ID = "Cray"
#==============================================================================

message(STATUS "-- Configuring Cray Fortran (ftn) compiler flags")

# Basic compiler flags
# -e m : Generate .mod files for modules
# -J   : Specify module output directory
set(FORTRAN_FLAGS_BASE "-e m -J ${CMAKE_Fortran_MODULE_DIRECTORY}")

# Add platform-specific preprocessor flags
set(FORTRAN_FLAGS_BASE "${FORTRAN_FLAGS_BASE} -D__linux__ -Dlinux")

# Debug configuration
# -e D : Enable all debugging features
# -R bcs : Enable bounds, conformance, and string checking
set(CMAKE_Fortran_FLAGS_DEBUG "-g -O0 -e D -R bcs ${FORTRAN_FLAGS_BASE}")

# Release flags - high optimization with vectorization
# -h fp3 : Aggressive floating-point optimizations (default)
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -h fp3 ${FORTRAN_FLAGS_BASE}")

# Reproducibility mode - strict IEEE 754 compliance
if(ENABLE_REPRODUCIBLE)
    # -h fp0 : Strict IEEE 754 compliance, no aggressive FP optimizations
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -h fp0 ${FORTRAN_FLAGS_BASE}")
    message(STATUS "-- Reproducibility: ENABLED (Cray -h fp0)")
endif()

# OpenMP parallelization - Cray uses -h omp
if(ENABLE_OPENMP)
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -h omp")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -h omp")
else()
    # Explicitly disable OpenMP if not requested
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -h noomp")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -h noomp")
endif()

# Report Cray CPU target if set
if(DEFINED ENV{CRAY_CPU_TARGET})
    message(STATUS "-- Cray CPU target: $ENV{CRAY_CPU_TARGET}")
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
