cmake_minimum_required(VERSION 3.16)

#==============================================================================
# Pre-Project Compiler Detection (prevents CMake cache loops)
#==============================================================================
if(NOT CMAKE_Fortran_COMPILER AND NOT DEFINED ENV{FC})
    # Only run detection if no compiler specified
    find_program(DEFAULT_FC NAMES gfortran f95 f77 fc fortran)
    if(DEFAULT_FC)
        get_filename_component(COMPILER_NAME ${DEFAULT_FC} NAME)
        get_filename_component(COMPILER_DIR ${DEFAULT_FC} DIRECTORY)
        
        # If we found a symlink (f95, f77, fc), try to find gfortran in same directory
        if(COMPILER_NAME MATCHES "^(f95|f77|fc|fortran)$")
            find_program(GFORTRAN_PREFERRED gfortran 
                PATHS ${COMPILER_DIR} 
                NO_DEFAULT_PATH
            )
            if(GFORTRAN_PREFERRED)
                set(CMAKE_Fortran_COMPILER ${GFORTRAN_PREFERRED} CACHE FILEPATH "Fortran compiler")
                message(STATUS "-- Pre-selected gfortran over ${COMPILER_NAME}: ${GFORTRAN_PREFERRED}")
            else()
                set(CMAKE_Fortran_COMPILER ${DEFAULT_FC} CACHE FILEPATH "Fortran compiler")
            endif()
            unset(GFORTRAN_PREFERRED CACHE)
        else()
            set(CMAKE_Fortran_COMPILER ${DEFAULT_FC} CACHE FILEPATH "Fortran compiler")
        endif()
        unset(DEFAULT_FC CACHE)
    endif()
endif()

#==============================================================================
# Read version from VERSION.md (single source of truth)
#==============================================================================
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION.md" _version_file)
string(REGEX MATCH "Version[^|]*\\|[ \t]*\\*\\*([0-9]+\\.[0-9]+\\.[0-9]+)\\*\\*" _ "${_version_file}")
if(NOT CMAKE_MATCH_1)
    message(FATAL_ERROR "Could not parse version from VERSION.md")
endif()
set(IOSEE_VERSION "${CMAKE_MATCH_1}")
message(STATUS "IOSEE version: ${IOSEE_VERSION} (from VERSION.md)")

project(iosee VERSION ${IOSEE_VERSION} LANGUAGES Fortran)

#==============================================================================
# Project Configuration
#==============================================================================
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(PROJECT_VERSION "${IOSEE_VERSION}")

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

#==============================================================================
# User-Configurable Options
#==============================================================================
option(ENABLE_OPENMP "Enable OpenMP parallelization" ON)
option(ENABLE_LTO "Enable Link-Time Optimization for Release builds" OFF)
option(ENABLE_FAST_MATH "Enable aggressive floating-point optimizations (may reduce accuracy)" OFF)
option(ENABLE_REPRODUCIBLE "Enable bit-for-bit reproducibility (strict IEEE 754)" OFF)
option(BUILD_TESTING "Build the test suite" OFF)

#==============================================================================
# Include Modules
#==============================================================================
include(cmake/SystemDetection.cmake)
include(cmake/BuildMessages.cmake)

# Print initial build message
print_build_header()

#==============================================================================
# Dependencies Detection
#==============================================================================
# OpenMP (conditional on ENABLE_OPENMP option)
if(ENABLE_OPENMP)
    find_package(OpenMP REQUIRED COMPONENTS Fortran)
    print_openmp_config()
    if(OpenMP_Fortran_FOUND)
        set(OMP_FLAGS ${OpenMP_Fortran_FLAGS})
        set(OMP_LIBS ${OpenMP_Fortran_LIBRARIES})

        if(NOT DEFINED ENV{OMP_NUM_THREADS})
            set(ENV{OMP_NUM_THREADS} ${DETECTED_CORES})
        endif()
    else()
        set(OMP_FLAGS "")
        set(OMP_LIBS "")
    endif()
else()
    message(STATUS "-- OpenMP: Disabled by user (ENABLE_OPENMP=OFF)")
    set(OMP_FLAGS "")
    set(OMP_LIBS "")
endif()

# NetCDF
include(cmake/FindNetCDF.cmake)
print_netcdf_config()

#==============================================================================
# Compiler Detection and Configuration
#==============================================================================

# Compiler information (detection already completed before project initialization)
get_filename_component(COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND NOT COMPILER_NAME STREQUAL "gfortran")
    message(STATUS "-- Using GNU Fortran compiler: ${CMAKE_Fortran_COMPILER}")
endif()

# Apply compiler-specific flags
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    include(cmake/CompilerFlagsGNU.cmake)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    include(cmake/CompilerFlagsIntel.cmake)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM")
    include(cmake/CompilerFlagsIntelLLVM.cmake)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    include(cmake/CompilerFlagsCray.cmake)
else()
    message(WARNING "Untested compiler: ${CMAKE_Fortran_COMPILER_ID}")
    set(CMAKE_Fortran_FLAGS "-O2")
    set(CMAKE_EXE_LINKER_FLAGS "")
endif()

# Print detailed configuration
print_configuration_summary()

#==============================================================================
# Build Directory Setup
#==============================================================================
# Note: Do NOT override CMAKE_BINARY_DIR - it prevents out-of-source builds
# and is managed by CMake automatically
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY})

#==============================================================================
# Source Organization
#==============================================================================
# Library sources (core modules)
set(CORE_SOURCES
    ${CMAKE_SOURCE_DIR}/lib/utils.F90
    ${CMAKE_SOURCE_DIR}/lib/configuration.F90
    ${CMAKE_SOURCE_DIR}/lib/netcdf_handler.F90
    ${CMAKE_SOURCE_DIR}/lib/mathlib.F90
    ${CMAKE_SOURCE_DIR}/lib/save_output.F90
)

# Main computation sources
set(SOLVER_SOURCES
    ${CMAKE_SOURCE_DIR}/src/sea_water_nk.F90
    ${CMAKE_SOURCE_DIR}/src/ocean_emissivity.F90
)

# Driver source
set(DRIVER_SOURCES
    ${CMAKE_SOURCE_DIR}/src/driver.F90
)

#==============================================================================
# Target Creation with Progress Messages
#==============================================================================
# Core library
add_library(ose_core ${CORE_SOURCES})
set_target_properties(ose_core PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Solver library  
add_library(ose_solver ${SOLVER_SOURCES})
set_target_properties(ose_solver PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)
target_link_libraries(ose_solver ose_core ${NETCDF_LIBS})

# Main executable
add_executable(iosee ${DRIVER_SOURCES})
set_target_properties(iosee PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Link all libraries
target_link_libraries(iosee ose_solver ose_core ${NETCDF_LIBS} ${OMP_LIBS})

# Include directories
target_include_directories(ose_core PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(ose_solver PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
target_include_directories(iosee PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Dependencies
add_dependencies(ose_solver ose_core)
add_dependencies(iosee ose_solver)

#==============================================================================
# Custom Targets
#==============================================================================
# Portable clean_all target (works on Windows, macOS, and Linux)
add_custom_target(clean_all
    COMMAND ${CMAKE_COMMAND} -E echo "Cleaning build files..."
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/build
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_SOURCE_DIR}/examples/iosee
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_SOURCE_DIR}/examples/data
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_SOURCE_DIR}/libose_core.a
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_SOURCE_DIR}/libose_solver.a
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_SOURCE_DIR}/cmake_install.cmake
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_SOURCE_DIR}/CMakeCache.txt
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_SOURCE_DIR}/Makefile
    COMMAND ${CMAKE_COMMAND} -E echo "Clean complete"
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Removing all build artifacts"
)

#==============================================================================
# Post-Build Actions
#==============================================================================
# Create examples directory and symbolic links
add_custom_command(TARGET iosee POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/examples
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            ../build/iosee
            ${CMAKE_SOURCE_DIR}/examples/iosee
    COMMAND ${CMAKE_COMMAND} -E create_symlink
            ../data
            ${CMAKE_SOURCE_DIR}/examples/data
    COMMENT "Creating symbolic links to executable and data in examples directory"
)

#==============================================================================
# Installation Rules
#==============================================================================
include(GNUInstallDirs)

# Install executable
install(TARGETS iosee
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install libraries
install(TARGETS ose_core ose_solver
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# Install Fortran modules
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/iosee
    FILES_MATCHING PATTERN "*.mod"
)

# Install data files
install(DIRECTORY ${CMAKE_SOURCE_DIR}/data/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/iosee
    FILES_MATCHING
        PATTERN "*.nc"
        PATTERN "*.dat"
)

# Install example configuration files
install(DIRECTORY ${CMAKE_SOURCE_DIR}/examples/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/iosee/examples
    FILES_MATCHING PATTERN "*.config"
)