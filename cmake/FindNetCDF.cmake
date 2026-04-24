#==============================================================================
# NetCDF Detection and Configuration
#==============================================================================

# Initialize variables
set(NETCDF_FOUND FALSE)
set(NETCDF_INCLUDE "")
set(NETCDF_LIBS "")

# Compiler-aware NetCDF detection for HPC environments
# Priority: Find NetCDF compiled with same compiler toolchain

# Detect current compiler toolchain
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(COMPILER_FAMILY "GCC")
    set(PREFERRED_NETCDF_PATTERNS "gcc" "gnu" "GCC" "GNU")
    set(AVOID_NETCDF_PATTERNS "intel" "icc" "ifort" "Intel" "iimpi")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(COMPILER_FAMILY "Intel")
    set(PREFERRED_NETCDF_PATTERNS "intel" "icc" "ifort" "Intel" "iimpi")
    set(AVOID_NETCDF_PATTERNS "gcc" "gnu" "GCC" "GNU")
else()
    set(COMPILER_FAMILY "Other")
    set(PREFERRED_NETCDF_PATTERNS "")
    set(AVOID_NETCDF_PATTERNS "")
endif()

# Try to find nf-config first (preferred method)
# Search for compiler-compatible nf-config first
if(PREFERRED_NETCDF_PATTERNS)
    # Look for nf-config in paths matching our compiler
    # Note: find_program() does not support PATH_REGEX, so we find all candidates
    # and then filter by pattern matching
    find_program(NF_CONFIG_CANDIDATE nf-config
        PATH_SUFFIXES bin
        PATHS
            "$ENV{NETCDF_ROOT}"
            "$ENV{NETCDF_HOME}"
            "$ENV{NETCDF_DIR}"
            "$ENV{NETCDF_FORTRAN_ROOT}"
            "/sw/eb/sw"
            "/opt"
            "/usr/local"
            "/apps"
        HINTS
            "$ENV{NETCDF_ROOT}/bin"
            "$ENV{NETCDF_HOME}/bin"
            "$ENV{NETCDF_DIR}/bin"
        NO_DEFAULT_PATH
    )

    # Check if found candidate matches preferred compiler patterns
    if(NF_CONFIG_CANDIDATE)
        foreach(pattern ${PREFERRED_NETCDF_PATTERNS})
            if(NF_CONFIG_CANDIDATE MATCHES "${pattern}")
                set(NF_CONFIG ${NF_CONFIG_CANDIDATE})
                message(STATUS "-- Found compiler-compatible nf-config: ${NF_CONFIG}")
                break()
            endif()
        endforeach()
    endif()
    unset(NF_CONFIG_CANDIDATE CACHE)
endif()

# Fallback: find any nf-config if no compiler-specific one found
if(NOT NF_CONFIG)
    find_program(NF_CONFIG nf-config)
    if(NF_CONFIG)
        # Check if this nf-config might be from wrong compiler
        foreach(avoid_pattern ${AVOID_NETCDF_PATTERNS})
            if(NF_CONFIG MATCHES "${avoid_pattern}")
                message(WARNING "-- Found NetCDF compiled with different compiler (${avoid_pattern}): ${NF_CONFIG}")
                message(WARNING "-- This may cause module compatibility issues with ${COMPILER_FAMILY} compiler")
                break()
            endif()
        endforeach()
    endif()
endif()

find_program(NC_CONFIG nc-config)

if(NF_CONFIG)
    execute_process(
        COMMAND ${NF_CONFIG} --fflags 
        OUTPUT_VARIABLE NETCDF_INCLUDE 
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND ${NF_CONFIG} --flibs 
        OUTPUT_VARIABLE NETCDF_LIBS 
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(NETCDF_FOUND TRUE)
    
elseif(NC_CONFIG)
    execute_process(
        COMMAND ${NC_CONFIG} --cflags 
        OUTPUT_VARIABLE NETCDF_INCLUDE 
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND ${NC_CONFIG} --libs 
        OUTPUT_VARIABLE NETCDF_LIBS 
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(NETCDF_FOUND TRUE)
    
else()
    # Fallback: Manual search for libraries
    
    if(APPLE)
        # Search Homebrew paths on macOS
        find_path(NETCDF_INCLUDE_DIR netcdf.mod
            PATHS 
                /opt/homebrew/include 
                /usr/local/include
                /opt/local/include
            PATH_SUFFIXES netcdf
        )
        
        find_library(NETCDF_FORTRAN_LIBRARY netcdff
            PATHS 
                /opt/homebrew/lib 
                /usr/local/lib
                /opt/local/lib
        )
        
        find_library(NETCDF_C_LIBRARY netcdf
            PATHS 
                /opt/homebrew/lib 
                /usr/local/lib
                /opt/local/lib
        )
        
        if(NETCDF_INCLUDE_DIR AND NETCDF_FORTRAN_LIBRARY AND NETCDF_C_LIBRARY)
            set(NETCDF_INCLUDE "-I${NETCDF_INCLUDE_DIR}")
            set(NETCDF_LIBS "${NETCDF_FORTRAN_LIBRARY} ${NETCDF_C_LIBRARY}")
            set(NETCDF_FOUND TRUE)
        endif()
        
    else()
        # Universal Linux/Unix fallback with HPC environment support
        # Prioritize compiler-compatible installations
        set(NETCDF_SEARCH_PATHS)
        set(NETCDF_LIB_PATHS)
        
        # Add preferred paths first (matching compiler toolchain)
        if(PREFERRED_NETCDF_PATTERNS)
            foreach(pattern ${PREFERRED_NETCDF_PATTERNS})
                list(APPEND NETCDF_SEARCH_PATHS
                    "/sw/eb/sw/netCDF-Fortran/*${pattern}*/include"
                    "/sw/eb/sw/netcdf-fortran/*${pattern}*/include"
                    "/opt/netcdf/*${pattern}*/include"
                    "/apps/netcdf/*${pattern}*/include"
                    "/usr/local/netcdf/*${pattern}*/include"
                )
                list(APPEND NETCDF_LIB_PATHS
                    "/sw/eb/sw/netCDF-Fortran/*${pattern}*/lib"
                    "/sw/eb/sw/netcdf-fortran/*${pattern}*/lib"
                    "/opt/netcdf/*${pattern}*/lib"
                    "/apps/netcdf/*${pattern}*/lib"
                    "/usr/local/netcdf/*${pattern}*/lib"
                    "/sw/eb/sw/netCDF-Fortran/*${pattern}*/lib64"
                    "/sw/eb/sw/netcdf-fortran/*${pattern}*/lib64"
                )
            endforeach()
        endif()
        
        # Add standard paths
        list(APPEND NETCDF_SEARCH_PATHS
            /usr/include
            /usr/local/include
            /opt/netcdf/include
            # Standard HPC paths
            /opt/netcdf-fortran/include
            /usr/local/netcdf/include
            /apps/netcdf/include
            /sw/netcdf/include
            # Module environment paths
            "$ENV{NETCDF_ROOT}/include"
            "$ENV{NETCDF_HOME}/include"
            "$ENV{NETCDF_DIR}/include"
            "$ENV{NETCDF_FORTRAN_ROOT}/include"
            "$ENV{NETCDF_FORTRAN_HOME}/include"
        )
        
        list(APPEND NETCDF_LIB_PATHS
            /usr/lib
            /usr/local/lib
            /opt/netcdf/lib
            /usr/lib64
            /usr/local/lib64
            # Standard HPC paths
            /opt/netcdf-fortran/lib
            /usr/local/netcdf/lib
            /apps/netcdf/lib
            /sw/netcdf/lib
            # Module environment paths
            "$ENV{NETCDF_ROOT}/lib"
            "$ENV{NETCDF_HOME}/lib"
            "$ENV{NETCDF_DIR}/lib"
            "$ENV{NETCDF_FORTRAN_ROOT}/lib"
            "$ENV{NETCDF_FORTRAN_HOME}/lib"
            "$ENV{NETCDF_ROOT}/lib64"
            "$ENV{NETCDF_HOME}/lib64"
            "$ENV{NETCDF_DIR}/lib64"
        )
        
        # Remove empty paths that come from undefined environment variables
        list(REMOVE_ITEM NETCDF_SEARCH_PATHS "ENV{}/include")
        list(REMOVE_ITEM NETCDF_LIB_PATHS "ENV{}/lib")
        list(REMOVE_ITEM NETCDF_LIB_PATHS "ENV{}/lib64")
        
        find_path(NETCDF_INCLUDE_DIR netcdf.mod
            PATHS ${NETCDF_SEARCH_PATHS}
            PATH_SUFFIXES netcdf fortran
        )
        
        find_library(NETCDF_FORTRAN_LIBRARY netcdff
            PATHS ${NETCDF_LIB_PATHS}
        )
        
        find_library(NETCDF_C_LIBRARY netcdf
            PATHS ${NETCDF_LIB_PATHS}
        )
        
        if(NETCDF_INCLUDE_DIR AND NETCDF_FORTRAN_LIBRARY AND NETCDF_C_LIBRARY)
            set(NETCDF_INCLUDE "-I${NETCDF_INCLUDE_DIR}")
            set(NETCDF_LIBS "${NETCDF_FORTRAN_LIBRARY} ${NETCDF_C_LIBRARY}")
            set(NETCDF_FOUND TRUE)
        endif()
    endif()
endif()

# Report results and provide guidance
if(NOT NETCDF_FOUND)
    message(FATAL_ERROR "NetCDF not found!")
    if(APPLE)
        message(FATAL_ERROR "On macOS, try: brew install netcdf-fortran")
    elseif(HPC_ENVIRONMENT)
        message(FATAL_ERROR "HPC System: NetCDF-Fortran not found or incompatible.")
        message(FATAL_ERROR "Current compiler: ${COMPILER_FAMILY} (${CMAKE_Fortran_COMPILER})")
        message(FATAL_ERROR "Solutions:")
        if(COMPILER_FAMILY STREQUAL "GCC")
            message(FATAL_ERROR "  1. Load GCC-compatible NetCDF module: 'module load netCDF-Fortran/X.X.X-GCC-X.X.X'")
            message(FATAL_ERROR "  2. Or switch to Intel toolchain: 'module purge && module load Intel/2023b netCDF-Fortran/4.6.1-iimpi-2023b'")
        elseif(COMPILER_FAMILY STREQUAL "Intel")
            message(FATAL_ERROR "  1. Load Intel-compatible NetCDF module: 'module load netCDF-Fortran/X.X.X-iimpi-XXXX'")
            message(FATAL_ERROR "  2. Or switch to GCC toolchain: 'module purge && module load GCC/X.X.X netCDF-Fortran/X.X.X-GCC-X.X.X'")
        else()
            message(FATAL_ERROR "  1. Load NetCDF-Fortran module matching your compiler")
            message(FATAL_ERROR "  2. Check 'module avail netCDF-Fortran' for available versions")
        endif()
    else()
        message(FATAL_ERROR "On Linux, try: sudo apt-get install libnetcdff-dev")
    endif()
elseif(NF_CONFIG)
    # Validate compiler compatibility  
    foreach(avoid_pattern ${AVOID_NETCDF_PATTERNS})
        if(NF_CONFIG MATCHES "${avoid_pattern}")
            message(WARNING "===== COMPILER COMPATIBILITY WARNING =====")
            message(WARNING "NetCDF compiled with different compiler toolchain!")
            message(WARNING "Current compiler: ${COMPILER_FAMILY} (${CMAKE_Fortran_COMPILER})")
            message(WARNING "NetCDF path suggests: ${avoid_pattern} compiler")
            message(WARNING "This may cause module compatibility errors during compilation.")
            message(WARNING "Recommended: Load NetCDF module matching your compiler toolchain.")
            message(WARNING "==========================================")
            break()
        endif()
    endforeach()
endif()