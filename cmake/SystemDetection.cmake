#==============================================================================
# Advanced System Detection Module
# Provides detailed hardware and platform information similar to satrs_lw_v1.0
#==============================================================================

# Get basic system information
# Robust Basic System Information with Fallbacks
cmake_host_system_information(RESULT PLATFORM QUERY OS_NAME)

# Core detection with fallbacks
cmake_host_system_information(RESULT DETECTED_CORES QUERY NUMBER_OF_LOGICAL_CORES)
cmake_host_system_information(RESULT DETECTED_PHYSICAL_CORES QUERY NUMBER_OF_PHYSICAL_CORES)

# Fallback core detection if CMake queries fail
if(NOT DETECTED_CORES OR DETECTED_CORES EQUAL 0)
    if(APPLE)
        execute_process(
            COMMAND sysctl -n hw.ncpu
            OUTPUT_VARIABLE DETECTED_CORES
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
    else()
        execute_process(
            COMMAND nproc
            OUTPUT_VARIABLE DETECTED_CORES
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(NOT DETECTED_CORES OR DETECTED_CORES EQUAL 0)
            # Last resort: check /proc/cpuinfo
            if(EXISTS "/proc/cpuinfo")
                execute_process(
                    COMMAND grep -c "^processor" /proc/cpuinfo
                    OUTPUT_VARIABLE DETECTED_CORES
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET
                )
            endif()
        endif()
    endif()
    # Set reasonable default if all detection methods fail
    if(NOT DETECTED_CORES OR DETECTED_CORES EQUAL 0)
        set(DETECTED_CORES 4)
        message(WARNING "-- Could not detect CPU cores, using default: ${DETECTED_CORES}")
    endif()
endif()

# Architecture detection with fallbacks
execute_process(COMMAND uname -m OUTPUT_VARIABLE ARCH OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT ARCH)
    # Fallback architecture detection
    if(CMAKE_SYSTEM_PROCESSOR)
        set(ARCH ${CMAKE_SYSTEM_PROCESSOR})
    else()
        set(ARCH "unknown")
        message(WARNING "-- Could not detect system architecture")
    endif()
endif()

# Universal HPC Environment Detection
set(HPC_ENVIRONMENT FALSE)
set(HPC_TYPE "Unknown")

# Check for common HPC module systems (universal across HPC platforms)
if(DEFINED ENV{MODULESHOME} OR DEFINED ENV{LMOD_DIR} OR DEFINED ENV{LMOD_VERSION} OR DEFINED ENV{MODULESPATH})
    set(HPC_ENVIRONMENT TRUE)
    if(DEFINED ENV{LMOD_VERSION})
        set(HPC_TYPE "LMOD")
    elseif(DEFINED ENV{MODULESHOME})
        set(HPC_TYPE "Environment Modules")
    elseif(DEFINED ENV{MODULESPATH})
        set(HPC_TYPE "Traditional Modules")
    endif()
endif()

# Check for HPC job schedulers (universal indicators)
if(PLATFORM MATCHES "Linux")
    if(DEFINED ENV{SLURM_CLUSTER_NAME} OR DEFINED ENV{SLURM_JOB_ID} OR DEFINED ENV{SLURM_PROCID})
        set(HPC_ENVIRONMENT TRUE)
        set(HPC_TYPE "SLURM Cluster")
    elseif(DEFINED ENV{PBS_JOBID} OR DEFINED ENV{PBS_O_WORKDIR})
        set(HPC_ENVIRONMENT TRUE)
        set(HPC_TYPE "PBS/Torque Cluster")
    elseif(DEFINED ENV{LSB_JOBID} OR DEFINED ENV{LSF_BINDIR})
        set(HPC_ENVIRONMENT TRUE)
        set(HPC_TYPE "LSF Cluster")
    elseif(DEFINED ENV{SGE_ROOT} OR DEFINED ENV{JOB_ID})
        set(HPC_ENVIRONMENT TRUE)
        set(HPC_TYPE "SGE/UGE Cluster")
    endif()
endif()

# Check for specific HPC systems and centers
if(DEFINED ENV{TACC_SYSTEM} OR DEFINED ENV{STOCKYARD})
    set(HPC_ENVIRONMENT TRUE)
    set(HPC_TYPE "TACC (${HPC_TYPE})")
elseif(DEFINED ENV{NERSC_HOST})
    set(HPC_ENVIRONMENT TRUE)
    set(HPC_TYPE "NERSC")
elseif(DEFINED ENV{OLCF_MODULEPATH_ROOT})
    set(HPC_ENVIRONMENT TRUE)
    set(HPC_TYPE "OLCF")
elseif(DEFINED ENV{ALCF_TOOLS})
    set(HPC_ENVIRONMENT TRUE)
    set(HPC_TYPE "ALCF")
endif()

# Additional HPC indicators for Linux systems
if(PLATFORM MATCHES "Linux" AND NOT HPC_ENVIRONMENT)
    # Check for common HPC directory structures
    if(EXISTS "/opt/modules" OR EXISTS "/usr/share/modules" OR EXISTS "/etc/modulefiles")
        set(HPC_ENVIRONMENT TRUE)
        set(HPC_TYPE "Module System Detected")
    # Check for common HPC software stacks
    elseif(EXISTS "/opt/intel" OR EXISTS "/opt/gcc" OR EXISTS "/apps" OR EXISTS "/sw")
        set(HPC_ENVIRONMENT TRUE)
        set(HPC_TYPE "HPC Software Stack")
    # Check for login node patterns (many HPC systems use these patterns)
    elseif(EXISTS "/scratch" OR EXISTS "/lustre" OR EXISTS "/gpfs")
        execute_process(
            COMMAND hostname
            OUTPUT_VARIABLE HOSTNAME_RESULT
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(HOSTNAME_RESULT MATCHES "(login|head|submit|frontend|gateway)")
            set(HPC_ENVIRONMENT TRUE)
            set(HPC_TYPE "HPC Login Node")
        endif()
    endif()
endif()

# Enhanced Apple System Detection (M-series and Intel Mac compatibility)
if(APPLE)
    # Detect architecture and CPU type
    if(ARCH MATCHES "arm64|aarch64")
        set(APPLE_SILICON TRUE)
        set(APPLE_INTEL FALSE)
        
        # Enhanced Apple Silicon chip detection
        execute_process(
            COMMAND sysctl -n machdep.cpu.brand_string
            OUTPUT_VARIABLE CPU_BRAND
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        
        if(CPU_BRAND MATCHES "Apple")
            # Comprehensive Apple chip detection including future chips
            if(CPU_BRAND MATCHES "M[0-9]+ Ultra")
                string(REGEX MATCH "M[0-9]+ Ultra" CHIP_TYPE "${CPU_BRAND}")
            elseif(CPU_BRAND MATCHES "M[0-9]+ Max")
                string(REGEX MATCH "M[0-9]+ Max" CHIP_TYPE "${CPU_BRAND}")
            elseif(CPU_BRAND MATCHES "M[0-9]+ Pro")
                string(REGEX MATCH "M[0-9]+ Pro" CHIP_TYPE "${CPU_BRAND}")
            elseif(CPU_BRAND MATCHES "M[0-9]+")
                string(REGEX MATCH "M[0-9]+" CHIP_TYPE "${CPU_BRAND}")
            else()
                set(CHIP_TYPE "Apple Silicon")
            endif()
            set(DETECTED_CPU_TYPE "Apple ${CHIP_TYPE}")
        else()
            set(DETECTED_CPU_TYPE "Apple Silicon")
        endif()
        
    elseif(ARCH MATCHES "x86_64")
        set(APPLE_SILICON FALSE)
        set(APPLE_INTEL TRUE)
        
        # Intel Mac detection with detailed CPU info
        execute_process(
            COMMAND sysctl -n machdep.cpu.brand_string
            OUTPUT_VARIABLE CPU_BRAND
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        
        if(CPU_BRAND)
            # Clean up Intel CPU name for display
            string(REGEX REPLACE "\\(R\\)" "" CPU_BRAND_CLEAN "${CPU_BRAND}")
            string(REGEX REPLACE "\\(TM\\)" "" CPU_BRAND_CLEAN "${CPU_BRAND_CLEAN}")
            string(REGEX REPLACE "  +" " " CPU_BRAND_CLEAN "${CPU_BRAND_CLEAN}")
            string(STRIP "${CPU_BRAND_CLEAN}" DETECTED_CPU_TYPE)
        else()
            set(DETECTED_CPU_TYPE "Intel Mac")
        endif()
        
    else()
        # Future-proofing for unknown Apple architectures
        set(APPLE_SILICON FALSE)
        set(APPLE_INTEL FALSE)
        set(DETECTED_CPU_TYPE "Apple ${ARCH}")
    endif()
    
else()
    set(APPLE_SILICON FALSE)
    if(APPLE)
        set(DETECTED_CPU_TYPE "Intel Mac")
    else()
        # Use CMake-native approach for CPU detection on Linux/HPC systems
        if(EXISTS "/proc/cpuinfo")
            file(READ "/proc/cpuinfo" CPUINFO_CONTENT)
            # Try to extract model name from cpuinfo using regex
            string(REGEX MATCH "model name[ \t]*:[ \t]*([^\n\r]+)" CPU_MATCH "${CPUINFO_CONTENT}")
            if(CPU_MATCH)
                set(CPU_BRAND "${CMAKE_MATCH_1}")
                string(STRIP "${CPU_BRAND}" DETECTED_CPU_TYPE)
            else()
                # Fallback: try processor info for different formats
                string(REGEX MATCH "processor[ \t]*:[ \t]*([^\n\r]+)" PROC_MATCH "${CPUINFO_CONTENT}")
                if(PROC_MATCH)
                    set(DETECTED_CPU_TYPE "Linux CPU")
                else()
                    if(HPC_ENVIRONMENT)
                        if(ARCH MATCHES "x86_64")
                            set(DETECTED_CPU_TYPE "HPC x86_64 CPU")
                        else()
                            set(DETECTED_CPU_TYPE "HPC ${ARCH} CPU")
                        endif()
                    else()
                        set(DETECTED_CPU_TYPE "Unknown CPU")
                    endif()
                endif()
            endif()
        else()
            # /proc/cpuinfo not available - likely HPC environment
            if(HPC_ENVIRONMENT)
                if(ARCH MATCHES "x86_64")
                    set(DETECTED_CPU_TYPE "HPC x86_64 CPU")
                elseif(ARCH MATCHES "aarch64")
                    set(DETECTED_CPU_TYPE "HPC ARM64 CPU")
                else()
                    set(DETECTED_CPU_TYPE "HPC ${ARCH} CPU")
                endif()
            else()
                set(DETECTED_CPU_TYPE "Unknown CPU")
            endif()
        endif()
    endif()
endif()

# Get system memory
if(APPLE)
    execute_process(
        COMMAND sysctl -n hw.memsize
        OUTPUT_VARIABLE MEMORY_BYTES
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(MEMORY_BYTES)
        math(EXPR MEMORY_GB "${MEMORY_BYTES} / 1024 / 1024 / 1024")
        set(DETECTED_MEMORY "${MEMORY_GB} GB")
    else()
        set(DETECTED_MEMORY "Unknown")
    endif()
else()
    # Use CMake-native approach for robust memory detection on Linux/HPC systems
    if(EXISTS "/proc/meminfo")
        file(READ "/proc/meminfo" MEMINFO_CONTENT)
        string(REGEX MATCH "MemTotal:[ \t]+([0-9]+)[ \t]+kB" MEMINFO_MATCH "${MEMINFO_CONTENT}")
        if(MEMINFO_MATCH)
            set(MEMORY_KB "${CMAKE_MATCH_1}")
            math(EXPR MEMORY_GB "${MEMORY_KB} / 1024 / 1024")
            set(DETECTED_MEMORY "${MEMORY_GB} GB")
        else()
            # Fallback: try simpler regex patterns for different /proc/meminfo formats
            string(REGEX MATCH "MemTotal:[ \t]*([0-9]+)" MEMINFO_MATCH_SIMPLE "${MEMINFO_CONTENT}")
            if(MEMINFO_MATCH_SIMPLE)
                set(MEMORY_KB "${CMAKE_MATCH_1}")
                if(MEMORY_KB AND MEMORY_KB GREATER 0)
                    math(EXPR MEMORY_GB "${MEMORY_KB} / 1024 / 1024")
                    set(DETECTED_MEMORY "${MEMORY_GB} GB")
                else()
                    set(DETECTED_MEMORY "Unknown")
                endif()
            else()
                # Additional fallback: try to find any number after MemTotal
                string(REGEX MATCH "MemTotal[^0-9]*([0-9]+)" MEMINFO_MATCH_LOOSE "${MEMINFO_CONTENT}")
                if(MEMINFO_MATCH_LOOSE)
                    set(MEMORY_KB "${CMAKE_MATCH_1}")
                    if(MEMORY_KB AND MEMORY_KB GREATER 0)
                        math(EXPR MEMORY_GB "${MEMORY_KB} / 1024 / 1024")
                        set(DETECTED_MEMORY "${MEMORY_GB} GB")
                    else()
                        set(DETECTED_MEMORY "Unknown")
                    endif()
                else()
                    set(DETECTED_MEMORY "Unknown")
                endif()
            endif()
        endif()
    else()
        # /proc/meminfo not available - try alternative methods
        set(DETECTED_MEMORY "Unknown")
        if(HPC_ENVIRONMENT)
        endif()
        
        # Additional fallback attempts for Linux systems
        if(NOT APPLE AND NOT HPC_ENVIRONMENT)
            execute_process(
                COMMAND free -m
                OUTPUT_VARIABLE FREE_OUTPUT
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET
            )
            if(FREE_OUTPUT)
                string(REGEX MATCH "Mem:[ \t]+([0-9]+)" FREE_MATCH "${FREE_OUTPUT}")
                if(FREE_MATCH)
                    set(MEMORY_MB "${CMAKE_MATCH_1}")
                    if(MEMORY_MB AND MEMORY_MB GREATER 0)
                        math(EXPR MEMORY_GB "${MEMORY_MB} / 1024")
                        set(DETECTED_MEMORY "${MEMORY_GB} GB")
                    endif()
                endif()
            endif()
        endif()
    endif()
endif()

# Detect optimal parallel job count for builds
if(NOT DEFINED PARALLEL_JOBS)
    if(APPLE)
        execute_process(
            COMMAND sysctl -n hw.ncpu
            OUTPUT_VARIABLE PARALLEL_JOBS
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
    else()
        execute_process(
            COMMAND nproc
            OUTPUT_VARIABLE PARALLEL_JOBS
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
    endif()
    if(NOT PARALLEL_JOBS OR PARALLEL_JOBS EQUAL 0)
        set(PARALLEL_JOBS 4)
    endif()
    message(STATUS "-- Parallel build jobs: ${PARALLEL_JOBS}")
endif()

# Universal Cache Optimization Settings
if(APPLE_SILICON)
    set(CACHE_OPTIMIZED TRUE)
    # Apple Silicon has unified memory architecture with large caches
    set(L1_CACHE_SIZE 128)    # KB (typically 128KB per P-core, 64KB per E-core)
    set(L2_CACHE_SIZE 12288)  # KB (12MB typical for Apple Silicon)
elseif(APPLE_INTEL)
    set(CACHE_OPTIMIZED TRUE)
    # Intel Mac systems have traditional cache hierarchy
    set(L1_CACHE_SIZE 64)     # KB (typical Intel L1)
    set(L2_CACHE_SIZE 512)    # KB (typical Intel L2 per core)
elseif(HPC_ENVIRONMENT)
    set(CACHE_OPTIMIZED TRUE)
    # HPC systems typically have larger caches
    if(ARCH MATCHES "x86_64")
        # Typical HPC x86_64 node configurations (Intel Xeon, AMD EPYC)
        set(L1_CACHE_SIZE 64)    # KB (per core)
        set(L2_CACHE_SIZE 1024)  # KB (1MB typical for modern HPC processors)
    elseif(ARCH MATCHES "aarch64")
        # ARM-based HPC systems (e.g., Fugaku, Archer2)
        set(L1_CACHE_SIZE 64)    # KB
        set(L2_CACHE_SIZE 512)   # KB
    else()
        # Conservative defaults for other HPC architectures
        set(L1_CACHE_SIZE 64)    # KB
        set(L2_CACHE_SIZE 512)   # KB
    endif()
else()
    # Standard desktop/workstation Linux systems
    set(CACHE_OPTIMIZED FALSE)
    set(L1_CACHE_SIZE 64)   # KB (typical for x86_64)
    set(L2_CACHE_SIZE 512)  # KB (typical for consumer processors)
endif()
