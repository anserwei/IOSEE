#==============================================================================
# Build Messaging Utilities
# Professional messaging system matching satrs_lw_v1.0 style
#==============================================================================

# Function to print initial build message
function(print_build_header)
    message(STATUS "Building Ocean Emissivity (${CMAKE_BUILD_TYPE} mode)...")
endfunction()

# Function to print OpenMP configuration
function(print_openmp_config)
    if(OpenMP_Fortran_FOUND)
    else()
        message(WARNING "-- OpenMP not found - parallel performance will be limited")
    endif()
endfunction()

# Function to print NetCDF configuration  
function(print_netcdf_config)
    if(NETCDF_FOUND)
        if(NF_CONFIG)
        elseif(NETCDF_INCLUDE_DIR)
        endif()
    else()
        message(FATAL_ERROR "-- NetCDF not found!")
    endif()
endfunction()

# Function to print build completion message
function(print_build_completion)
    message(STATUS "Ocean Emissivity built with OpenMP support: ${OpenMP_Fortran_FOUND}")
    message(STATUS "OpenMP thread count: ${DETECTED_CORES}")
    message(STATUS "Run with ./run_optimized.sh for best performance")
    message(STATUS "Build complete!")
    message(STATUS "======================  How to run  =======================")
    message(STATUS "Manual configuration method:")
    message(STATUS "  ./ose_test config.file <output_file>")
    message(STATUS "=========================================================")
endfunction()

# Function to print detailed configuration summary
function(print_configuration_summary)
    message(STATUS "")
    if(HPC_ENVIRONMENT)
    endif()
    if(APPLE_SILICON)
    elseif(APPLE_INTEL)
    elseif(HPC_ENVIRONMENT)
    endif()
    if(DETECTED_PHYSICAL_CORES STREQUAL DETECTED_CORES)
        message(STATUS "-- Building with ${DETECTED_CORES} cores...")
    else()
        message(STATUS "-- Building with ${DETECTED_CORES} cores (${DETECTED_PHYSICAL_CORES} physical cores)...")
    endif()
    message(STATUS "")
endfunction()