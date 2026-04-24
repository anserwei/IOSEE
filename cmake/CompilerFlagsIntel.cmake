#==============================================================================
# Intel Fortran Compiler Configuration
#==============================================================================


# Basic compiler flags with platform-specific preprocessor definitions
set(BASIC_FLAGS "-fpp -assume realloc_lhs -standard-semantics")

# Add platform-specific preprocessor flags for proper platform detection
if(UNIX AND NOT APPLE)
    # Linux/Unix systems (including HPC)
    set(BASIC_FLAGS "${BASIC_FLAGS} -D__linux__ -Dlinux")
elseif(APPLE)
    # Apple systems
    set(BASIC_FLAGS "${BASIC_FLAGS} -D__APPLE__")
endif()

# Intel-specific aggressive optimizations for maximum performance
set(OPT_FLAGS "-O3 -xHost -mtune=native")
set(OPT_FLAGS "${OPT_FLAGS} -ip -ipo -no-prec-div -fp-model fast=2")
set(OPT_FLAGS "${OPT_FLAGS} -qopenmp-simd -qopt-prefetch=5")
set(OPT_FLAGS "${OPT_FLAGS} -qopt-streaming-stores always")
set(OPT_FLAGS "${OPT_FLAGS} -qopt-streaming-cache-evict=0")
set(OPT_FLAGS "${OPT_FLAGS} -qopt-multiple-gather-scatter-by-shuffles")
set(OPT_FLAGS "${OPT_FLAGS} -qopt-zmm-usage=high -qopt-dynamic-align")
set(OPT_FLAGS "${OPT_FLAGS} -parallel -par-threshold0 -qoverride-limits")

# Reproducibility mode - strict IEEE 754 compliance
if(ENABLE_REPRODUCIBLE)
    string(REPLACE "-fp-model fast=2" "-fp-model strict" OPT_FLAGS "${OPT_FLAGS}")
    string(REPLACE "-no-prec-div" "" OPT_FLAGS "${OPT_FLAGS}")
    set(OPT_FLAGS "${OPT_FLAGS} -fimf-arch-consistency=true")
    message(STATUS "-- Reproducibility: ENABLED (strict IEEE 754)")
endif()

# Warning flags
set(WARN_FLAGS "-warn all")

# Debug configuration
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(DEBUG_FLAGS "-check bounds -traceback -g -O0")
    set(OPT_FLAGS "-O1")
    set(LTO_FLAGS "")
else()
    set(DEBUG_FLAGS "")
    set(LTO_FLAGS "-ipo -qipo")
endif()

# Combine all flags
set(CMAKE_Fortran_FLAGS "${BASIC_FLAGS} ${OPT_FLAGS} ${WARN_FLAGS} ${DEBUG_FLAGS} ${OMP_FLAGS} ${NETCDF_INCLUDE}")

# Linker flags
if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LTO_FLAGS}")
endif()
