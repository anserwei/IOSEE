! configuration.F90
!
! Purpose:
!   Configuration management and automatic hardware detection for IOSEE.
!   Provides Fortran namelist parsing, input validation, and hardware-aware
!   performance parameter optimization.
!
! This module handles all user-configurable parameters for ocean emissivity
! calculations, including spectral ranges, environmental conditions, and
! computational settings.
!
! Features:
!   - Fortran namelist parsing with comprehensive validation
!   - Automatic CPU architecture detection (Apple Silicon, Intel/AMD, ARM)
!   - OpenMP thread detection with HPC job scheduler awareness
!   - Performance parameter auto-tuning based on problem size
!   - Array capacity management (up to 200 elements per parameter)
!   - Multi-mode support: spectral, narrow-band, polarized
!   - Output control: no_multiple_reflection option
!
! Configuration Parameters:
!   - wavenum_start, wavenum_end, n_wavenum: Spectral range
!   - sea_surface_temperature: SST in Kelvin
!   - sea_surface_salinity: Salinity in PSU
!   - view_angle: Viewing angles from nadir (degrees)
!   - wind: Wind speeds at 10m height (m/s)
!   - mode: 'spectral' or 'narrow-band'
!   - polarization: 'enable' or 'disable'
!   - no_multiple_reflection: Output first-order emissivity only
!
! Public Interface:
!   - config_type: Main configuration derived type
!     - initialize(): Hardware detection and auto-tuning
!     - read(): Namelist file parsing with validation
!     - validate(): Input parameter validation
!     - print_summary(): Configuration summary output
!
! (C) Copyright 2025-, Texas A&M University.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! Maintainer: Atmospheric & Oceanic Optics Group, Department of Atmospheric Sciences, Texas A&M University
! Author:     Dr. Jian Wei
! Email:      anser@tamu.edu
!
!
module configuration
    use utils, only : jprd, jpim
    use omp_lib, only : omp_get_max_threads
    implicit none
    public

    !---------------------------------------------------------------------
    ! Derived type containing all the configuration information.
    !---------------------------------------------------------------------
    type config_type
        ! user-settable parameters

        !--------------------------
        ! Input data file for water refractive index
        !--------------------------
        character(len=511) :: refractive_index_filename = "data/water_optical_constants.nc"
        
        !--------------------------
        ! Calculation parameters
        !--------------------------
        ! Mode: 'spectral' or 'narrow-band'
        character(len=20)  :: mode = 'spectral'
        
        ! Arrays for multiple bands (up to 851 bands)
        real(kind=jprd)    :: wavenum_start(851) = 0.0_jprd   ! Starting wavenumbers (cm-1)
        real(kind=jprd)    :: wavenum_end(851)   = 0.0_jprd   ! Ending wavenumbers (cm-1)
        integer(kind=jpim) :: n_wavenum(851)     = 0          ! Number of wavenumbers per band
        integer(kind=jpim) :: n_bands            = 1          ! Number of bands (derived from input)
        
        ! Wavelength parameters (alternative to wavenumber)
        real(kind=jprd)    :: wavelen_min(851)   = 0.0_jprd   ! Starting wavelengths (um)
        real(kind=jprd)    :: wavelen_max(851)   = 0.0_jprd   ! Ending wavelengths (um)
        integer(kind=jpim) :: n_wavelen(851)     = 0          ! Number of wavelengths per band
        logical            :: is_wavelength_input = .false.    ! Flag to indicate wavelength input
        logical            :: no_multiple_reflection = .false. ! Output emissivity without multiple reflection
        logical            :: use_effective_angle = .true.    ! Use effective view angle LUT

        ! Environmental parameters
        real(kind=jprd)    :: sea_surface_temperature   = 298.0_jprd     ! Sea surface temperature (K)
        real(kind=jprd)    :: sea_surface_salinity      = 35.0_jprd      ! Sea surface salinity (PSU)
        
        ! Viewing angles and wind speeds
        real(kind=jprd)    :: view_angle(851)  = 0.0_jprd     ! Viewing angles (degrees)
        real(kind=jprd)    :: winds(851)       = 0.0_jprd     ! Wind speeds (m/s)
        integer(kind=jpim) :: n_angles        = 0            ! Number of viewing angles (derived)
        integer(kind=jpim) :: n_winds         = 0            ! Number of wind speeds (derived)
        
        ! Polarization parameters
        character(len=20)  :: polarization    = 'disable'     ! 'disable' or 'enable'

        ! Output file configuration
        character(len=512) :: output_file     = ""            ! User-definable output filename (empty = use default)

        !--------------------------
        ! Hardware detection and performance parameters
        !--------------------------
        ! Automatically detected hardware parameters
        integer(kind=jpim) :: MAX_THREADS          = 1      ! Auto-detected
        logical            :: is_apple_silicon     = .false. ! Auto-detected
        character(len=64)  :: cpu_architecture     = "unknown"
        
        ! Cache and vectorization parameters (auto-tuned based on hardware)
        integer(kind=jpim) :: CACHE_LINE_SIZE      = 64      ! Will be auto-detected
        integer(kind=jpim) :: VECTOR_LENGTH        = 8       ! Will be auto-detected
        integer(kind=jpim) :: L1_CACHE_SIZE        = 32768   ! Will be auto-detected
        integer(kind=jpim) :: L2_BLOCK_SIZE        = 256     ! Will be auto-tuned
        integer(kind=jpim) :: MIN_CHUNK_SIZE       = 32      ! Will be auto-tuned
        integer(kind=jpim) :: MAX_CHUNK_SIZE       = 512     ! Will be auto-tuned
        integer(kind=jpim) :: MIN_PARALLEL_SIZE    = 1000    ! Will be auto-tuned
        integer(kind=jpim) :: MIN_WORK_PER_THREAD  = 64      ! Will be auto-tuned
        
        ! Performance optimization flags
        logical            :: use_guided_scheduling = .true.
        logical            :: use_simd_optimization = .true.
        logical            :: use_cache_blocking    = .true.

    contains
        procedure :: read => read_configuration
        procedure :: initialize => initialize_hardware_config
        procedure :: print_config => print_configuration
        procedure :: optimize_for_hardware => optimize_for_hardware
    end type config_type

contains

    !-----------------------------------------------------------------------
    ! Subroutine: read_configuration
    !
    ! Reads user configuration from a namelist file or unit. Overwrites
    ! the default values if present in the file. 
    !-----------------------------------------------------------------------
    subroutine read_configuration(this, file_name, unit, is_success)
        implicit none
        class(config_type), intent(inout)         :: this
        character(len=511), intent(in),  optional :: file_name
        integer,           intent(in),  optional :: unit
        logical,           intent(out), optional :: is_success

        integer            :: iosopen, iosread
        integer            :: iunit
        integer(kind=jpim) :: i  ! Loop variable for validation

        ! Local copies for namelist read
        character(len=511) :: auxiliary_filename          ! Old name (backward compat)
        character(len=511) :: refractive_index_filename   ! New preferred name
        integer(kind=jpim) :: CACHE_LINE_SIZE
        integer(kind=jpim) :: VECTOR_LENGTH
        integer(kind=jpim) :: L1_CACHE_SIZE
        integer(kind=jpim) :: L2_BLOCK_SIZE
        integer(kind=jpim) :: MIN_CHUNK_SIZE
        integer(kind=jpim) :: MAX_CHUNK_SIZE
        integer(kind=jpim) :: MIN_PARALLEL_SIZE
        integer(kind=jpim) :: MIN_WORK_PER_THREAD
        
        ! Calculation parameters - support both scalar and array formats
        character(len=20)  :: mode
        character(len=20)  :: polarization           ! Polarization parameter
        real(kind=jprd)    :: wavenum_start(851)
        real(kind=jprd)    :: wavenum_end(851)
        integer(kind=jpim) :: n_wavenum(851)
        real(kind=jprd)    :: temperature             ! For backward compatibility
        real(kind=jprd)    :: salinity                ! For backward compatibility
        real(kind=jprd)    :: sea_surface_temperature ! New name
        real(kind=jprd)    :: sea_surface_salinity    ! New name
        real(kind=jprd)    :: view_angle(851)
        real(kind=jprd)    :: winds(851)
        real(kind=jprd)    :: wind_values(851)  ! Alternative name for winds
        real(kind=jprd)    :: wind(851)         ! New consistent name
        
        ! Range-based parameters (Option 2)
        real(kind=jprd)    :: view_angle_start
        real(kind=jprd)    :: view_angle_end
        integer(kind=jpim) :: n_view_angle
        real(kind=jprd)    :: wind_start
        real(kind=jprd)    :: wind_end
        integer(kind=jpim) :: n_wind
        
        ! Scalar versions for backward compatibility
        real(kind=jprd)    :: wavenum_start_scalar
        real(kind=jprd)    :: wavenum_end_scalar
        integer(kind=jpim) :: n_wavenum_scalar
        
        ! Wavelength parameters
        real(kind=jprd)    :: wavelen_min(851)
        real(kind=jprd)    :: wavelen_max(851)
        integer(kind=jpim) :: n_wavelen(851)

        ! Output control
        logical            :: no_multiple_reflection
        logical            :: use_effective_angle
        character(len=512) :: output_file

        namelist /configuration/ auxiliary_filename, refractive_index_filename, &
            CACHE_LINE_SIZE, &
            VECTOR_LENGTH, &
            L1_CACHE_SIZE, &
            L2_BLOCK_SIZE, &
            MIN_CHUNK_SIZE, &
            MAX_CHUNK_SIZE, &
            MIN_PARALLEL_SIZE, &
            MIN_WORK_PER_THREAD, &
            mode, &
            polarization, &
            wavenum_start, &
            wavenum_end, &
            n_wavenum, &
            wavenum_start_scalar, &
            wavenum_end_scalar, &
            n_wavenum_scalar, &
            temperature, &
            salinity, &
            sea_surface_temperature, &
            sea_surface_salinity, &
            view_angle, &
            view_angle_start, &
            view_angle_end, &
            n_view_angle, &
            winds, &
            wind_values, &
            wind, &
            wind_start, &
            wind_end, &
            n_wind, &
            wavelen_min, &
            wavelen_max, &
            n_wavelen, &
            no_multiple_reflection, &
            use_effective_angle, &
            output_file

        if (.not. present(file_name) .and. .not. present(unit)) then
            write(*,'(a)') "===== Error: missing configuration file argument!"
            if (present(is_success)) is_success = .false.
            return
        end if
        
        ! Initialize local variables with defaults
        auxiliary_filename = this%refractive_index_filename
        refractive_index_filename = this%refractive_index_filename
        CACHE_LINE_SIZE = this%CACHE_LINE_SIZE
        VECTOR_LENGTH = this%VECTOR_LENGTH
        L1_CACHE_SIZE = this%L1_CACHE_SIZE
        L2_BLOCK_SIZE = this%L2_BLOCK_SIZE
        MIN_CHUNK_SIZE = this%MIN_CHUNK_SIZE
        MAX_CHUNK_SIZE = this%MAX_CHUNK_SIZE
        MIN_PARALLEL_SIZE = this%MIN_PARALLEL_SIZE
        MIN_WORK_PER_THREAD = this%MIN_WORK_PER_THREAD
        mode = this%mode
        polarization = this%polarization
        wavenum_start = this%wavenum_start
        wavenum_end = this%wavenum_end
        n_wavenum = this%n_wavenum
        ! Initialize scalar defaults
        wavenum_start_scalar = this%wavenum_start(1)
        wavenum_end_scalar = this%wavenum_end(1)
        n_wavenum_scalar = this%n_wavenum(1)
        temperature = this%sea_surface_temperature  ! Support old name
        salinity = this%sea_surface_salinity        ! Support old name
        sea_surface_temperature = this%sea_surface_temperature
        sea_surface_salinity = this%sea_surface_salinity
        view_angle = this%view_angle
        winds = this%winds
        wind_values = this%winds  ! Initialize both with same default
        wind = this%winds         ! Initialize with same default
        
        ! Initialize wavelength parameters with default values
        wavelen_min = this%wavelen_min
        wavelen_max = this%wavelen_max
        n_wavelen = this%n_wavelen

        ! Initialize output control parameters
        no_multiple_reflection = this%no_multiple_reflection
        use_effective_angle = this%use_effective_angle
        output_file = this%output_file

        ! Initialize range-based parameters with default values
        view_angle_start = 0.0_jprd
        view_angle_end = 0.0_jprd
        n_view_angle = 0
        wind_start = 0.0_jprd
        wind_end = 0.0_jprd
        n_wind = 0

        if (present(file_name)) then
            iunit = 8
            open(unit=iunit, file=trim(file_name), iostat=iosopen, recl=512)
        else
            iosopen = 0
            iunit   = unit
        end if

        if (iosopen == 0) then
            read(unit=iunit, nml=configuration, iostat=iosread)
            if (iosread /= 0) then
                if (present(is_success)) then
                    is_success = .false.
                else if (present(file_name)) then
                    write(*,'(a,a,a)') "===== Error reading namelist from file: '", &
                                       trim(file_name), "'"
                    write(*,'(a,i0)') "===== iostat=", iosread
                    if (iosread == 5010) then
                        write(*,'(a)') "===== Error 5010: Invalid namelist format detected"
                        write(*,'(a)') "===== Check for proper Fortran namelist syntax"
                    end if
                    close(unit=iunit)
                else
                    write(*,'(a,i0)') "===== Error reading namelist from unit ", iunit
                end if
            else
                write(*,'(a)') "Configuration file read successfully"
            end if
            if (present(file_name)) close(unit=iunit)
        else
            if (present(is_success)) is_success = .false.
            return
        end if

        if (present(is_success)) is_success = .true.

        ! Copy namelist values into 'this' config
        ! Handle both old and new variable names for backward compatibility
        if (refractive_index_filename /= this%refractive_index_filename) then
            this%refractive_index_filename = refractive_index_filename
        else if (auxiliary_filename /= this%refractive_index_filename) then
            this%refractive_index_filename = auxiliary_filename
        end if
        this%CACHE_LINE_SIZE         = CACHE_LINE_SIZE
        this%VECTOR_LENGTH           = VECTOR_LENGTH
        this%L1_CACHE_SIZE           = L1_CACHE_SIZE
        this%L2_BLOCK_SIZE           = L2_BLOCK_SIZE
        this%MIN_CHUNK_SIZE          = MIN_CHUNK_SIZE
        this%MAX_CHUNK_SIZE          = MAX_CHUNK_SIZE
        this%MIN_PARALLEL_SIZE       = MIN_PARALLEL_SIZE
        this%MIN_WORK_PER_THREAD     = MIN_WORK_PER_THREAD
        this%mode                    = mode
        this%polarization            = polarization
        ! Handle both old and new variable names for backward compatibility
        if (sea_surface_temperature /= this%sea_surface_temperature) then
            this%sea_surface_temperature = sea_surface_temperature
        else if (temperature /= this%sea_surface_temperature) then
            this%sea_surface_temperature = temperature
        end if
        
        if (sea_surface_salinity /= this%sea_surface_salinity) then
            this%sea_surface_salinity = sea_surface_salinity
        else if (salinity /= this%sea_surface_salinity) then
            this%sea_surface_salinity = salinity
        end if
        ! Handle view angle parameters - check for range-based format first, then arrays
        if (n_view_angle > 0 .and. view_angle_end > view_angle_start) then
            ! Range-based view angle specification (Option 2)
            this%n_angles = n_view_angle
            if (n_view_angle == 1) then
                this%view_angle(1) = view_angle_start
            else
                ! Generate linear spacing: view_angle = linspace(view_angle_start, view_angle_end, n_view_angle)
                do i = 1, n_view_angle
                    this%view_angle(i) = view_angle_start + (view_angle_end - view_angle_start) * &
                                         real(i - 1, jprd) / real(n_view_angle - 1, jprd)
                end do
            end if
            ! Zero out unused array elements
            this%view_angle(n_view_angle+1:) = 0.0_jprd
        else
            ! Array-based view angle specification (Option 1)
            this%view_angle = view_angle
        end if
        
        ! Handle scalar vs array format for wavenumber parameters
        if (wavenum_start_scalar /= 0.0_jprd) then
            ! Scalar format detected, convert to array format
            this%wavenum_start(1) = wavenum_start_scalar
            this%wavenum_end(1)   = wavenum_end_scalar
            this%n_wavenum(1)     = n_wavenum_scalar
            this%wavenum_start(2:) = 0.0_jprd
            this%wavenum_end(2:)   = 0.0_jprd
            this%n_wavenum(2:)     = 0
        else
            ! Array format
            this%wavenum_start = wavenum_start
            this%wavenum_end   = wavenum_end
            this%n_wavenum     = n_wavenum
        end if
        
        ! Copy output control parameters
        this%no_multiple_reflection = no_multiple_reflection
        this%use_effective_angle = use_effective_angle
        this%output_file = output_file

        ! Handle wavelength parameters and convert to wavenumber if needed
        ! Detect wavelength input by checking if any wavelength parameters are non-zero
        if (any(n_wavelen > 0) .or. any(wavelen_min > 0.0_jprd) .or. any(wavelen_max > 0.0_jprd)) then
            this%is_wavelength_input = .true.
            this%wavelen_min = wavelen_min
            this%wavelen_max = wavelen_max
            this%n_wavelen = n_wavelen
            
            ! Convert wavelength to wavenumber for internal processing
            ! wavenumber (cm^-1) = 10000 / wavelength (um)
            this%n_bands = count(this%n_wavelen > 0)
            
            do i = 1, this%n_bands
                if (this%n_wavelen(i) > 0) then
                    ! Convert wavelength bounds to wavenumber bounds
                    ! Note: wavelength min -> wavenumber max and vice versa (inverse relationship)
                    this%wavenum_start(i) = 10000.0_jprd / this%wavelen_max(i)
                    this%wavenum_end(i)   = 10000.0_jprd / this%wavelen_min(i)
                    this%n_wavenum(i)     = this%n_wavelen(i)
                end if
            end do
            
            ! Zero out unused array elements
            if (this%n_bands < 851) then
                this%wavenum_start(this%n_bands+1:) = 0.0_jprd
                this%wavenum_end(this%n_bands+1:)   = 0.0_jprd
                this%n_wavenum(this%n_bands+1:)     = 0
                this%wavelen_min(this%n_bands+1:)   = 0.0_jprd
                this%wavelen_max(this%n_bands+1:)   = 0.0_jprd
                this%n_wavelen(this%n_bands+1:)     = 0
            end if
        else
            this%is_wavelength_input = .false.
            ! Zero out wavelength arrays for wavenumber input
            this%wavelen_min = 0.0_jprd
            this%wavelen_max = 0.0_jprd
            this%n_wavelen   = 0
        end if
        
        ! Handle wind parameters - check for range-based format first, then arrays
        if (n_wind > 0 .and. wind_end > wind_start) then
            ! Range-based wind specification (Option 2)
            this%n_winds = n_wind
            if (n_wind == 1) then
                this%winds(1) = wind_start
            else
                ! Generate linear spacing: winds = linspace(wind_start, wind_end, n_wind)
                do i = 1, n_wind
                    this%winds(i) = wind_start + (wind_end - wind_start) * &
                                   real(i - 1, jprd) / real(n_wind - 1, jprd)
                end do
            end if
            ! Zero out unused array elements
            this%winds(n_wind+1:) = 0.0_jprd
        else
            ! Array-based wind specification (Option 1) - prioritize 'wind' for consistency
            if (any(wind /= 0.0_jprd)) then
                this%winds = wind
            else if (any(wind_values /= 0.0_jprd)) then
                this%winds = wind_values
            else
                this%winds = winds
            end if
        end if
        
        ! Determine number of bands, angles, and winds from arrays
        this%n_bands = count(this%n_wavenum > 0)
        
        ! For angles, count valid values (unless using range-based specification)
        if (n_view_angle > 0) then
            ! Already set above for range-based
        else
            ! For view angles, find the last consecutive valid entry from the start
            this%n_angles = 0
            do i = 1, 851
                ! Check if this is a valid view angle entry
                if (i == 1) then
                    ! First entry: accept any valid angle (including 0.0)
                    if (this%view_angle(i) >= 0.0_jprd .and. this%view_angle(i) <= 85.0_jprd) then
                        this%n_angles = 1
                    else
                        exit
                    end if
                else
                    ! Subsequent entries: continue if non-zero and within valid range
                    if (this%view_angle(i) > 0.0_jprd .and. this%view_angle(i) <= 85.0_jprd) then
                        this%n_angles = i
                    else
                        exit
                    end if
                end if
            end do
        end if
        
        ! For winds, count non-zero values, but if first value is 0, assume it is a valid wind speed
        if (n_wind > 0) then
            ! Already set above for range-based
        else
            if (this%winds(1) == 0.0_jprd) then
                this%n_winds = count(this%winds /= 0.0_jprd) + 1  ! Include the zero wind speed
            else
                this%n_winds = count(this%winds /= 0.0_jprd)
            end if
        end if
        
        ! Validate array size limits (maximum 851 elements)
        if (this%n_angles > 851) then
            write(*,'(a)') "======================================"
            write(*,'(a)') "ERROR: Array size limit exceeded!"
            write(*,'(a,i0,a)') "Number of view angles (", this%n_angles, ") exceeds maximum limit of 851"
            write(*,'(a)') "Please reduce the number of view angles in your configuration file"
            write(*,'(a)') "======================================"
            if (present(is_success)) then
                is_success = .false.
                return
            else
                stop 1
            end if
        end if
        
        if (this%n_winds > 851) then
            write(*,'(a)') "======================================"
            write(*,'(a)') "ERROR: Array size limit exceeded!"
            write(*,'(a,i0,a)') "Number of wind speeds (", this%n_winds, ") exceeds maximum limit of 851"
            write(*,'(a)') "Please reduce the number of wind speeds in your configuration file"
            write(*,'(a)') "======================================"
            if (present(is_success)) then
                is_success = .false.
                return
            else
                stop 1
            end if
        end if
        
        if (this%n_bands > 851) then
            write(*,'(a)') "======================================"
            write(*,'(a)') "ERROR: Array size limit exceeded!"
            write(*,'(a,i0,a)') "Number of bands (", this%n_bands, ") exceeds maximum limit of 851"
            write(*,'(a)') "Please reduce the number of bands in your configuration file"
            write(*,'(a)') "======================================"
            if (present(is_success)) then
                is_success = .false.
                return
            else
                stop 1
            end if
        end if
        
        ! Validate wind speeds: must be in range [0, 18] m/s
        do i = 1, this%n_winds
            if (this%winds(i) < 0.0_jprd .or. this%winds(i) > 18.0_jprd) then
                write(*,'(a,f6.2,a)') "ERROR: Wind speed ", this%winds(i), " m/s is out of valid range [0, 18] m/s"
                write(*,'(a)') "Valid wind speed range: 0.0 to 18.0 m/s"
                if (present(is_success)) then
                    is_success = .false.
                    return
                else
                    stop 1
                end if
            end if
        end do
        
        ! Validate sea surface temperature: must be in range [271, 311] K
        if (this%sea_surface_temperature < 271.0_jprd .or. this%sea_surface_temperature > 311.0_jprd) then
            write(*,'(a,f7.2,a)') "ERROR: Sea surface temperature ", this%sea_surface_temperature, &
                " K is out of valid range [271, 311] K"
            write(*,'(a)') "Valid sea surface temperature range: 271.0 to 311.0 K"
            if (present(is_success)) then
                is_success = .false.
                return
            else
                stop 1
            end if
        end if

        ! Validate view angles: must be in range [0, 85] degrees
        do i = 1, this%n_angles
            if (this%view_angle(i) < 0.0_jprd .or. this%view_angle(i) > 85.0_jprd) then
                write(*,'(a,f6.2,a)') "ERROR: View angle ", this%view_angle(i), " degrees is out of valid range [0, 85] degrees"
                write(*,'(a)') "Valid view angle range: 0.0 to 85.0 degrees"
                if (present(is_success)) then
                    is_success = .false.
                    return
                else
                    stop 1
                end if
            end if
        end do
        
        ! Validate wavenumber: must be in range [10, 5000] cm^-1
        ! This validation ensures that the wavenumber range is within the physically meaningful
        ! infrared spectrum and within the bounds of available water optical constants data.
        ! The range 10-5000 cm^-1 corresponds to wavelengths of 2-1000 μm.
        do i = 1, this%n_bands
            if (this%n_wavenum(i) > 0) then
                if (this%wavenum_start(i) < 10.0_jprd .or. this%wavenum_start(i) > 5000.0_jprd) then
                    write(*,'(a,f8.2,a)') "ERROR: Wavenumber start ", this%wavenum_start(i), &
                        " cm^-1 is out of valid range [10, 5000] cm^-1"
                    write(*,'(a)') "Valid wavenumber range: 10.0 to 5000.0 cm^-1"
                    if (present(is_success)) then
                        is_success = .false.
                        return
                    else
                        stop 1
                    end if
                end if
                if (this%wavenum_end(i) < 10.0_jprd .or. this%wavenum_end(i) > 5000.0_jprd) then
                    write(*,'(a,f8.2,a)') "ERROR: Wavenumber end ", this%wavenum_end(i), &
                        " cm^-1 is out of valid range [10, 5000] cm^-1"
                    write(*,'(a)') "Valid wavenumber range: 10.0 to 5000.0 cm^-1"
                    if (present(is_success)) then
                        is_success = .false.
                        return
                    else
                        stop 1
                    end if
                end if
            end if
        end do
        
        ! Validate wavelength: must be in range [2, 1000] μm
        ! This validation is only performed when wavelength input is detected.
        ! The range 2-1000 μm corresponds to the infrared spectrum used in ocean emissivity calculations.
        ! This range matches the wavenumber range of 10-5000 cm^-1 (since wavenumber = 10000/wavelength).
        if (this%is_wavelength_input) then
            do i = 1, this%n_bands
                if (this%n_wavelen(i) > 0) then
                    if (this%wavelen_min(i) < 2.0_jprd .or. this%wavelen_min(i) > 1000.0_jprd) then
                        write(*,'(a,f8.2,a)') "ERROR: Wavelength min ", this%wavelen_min(i), &
                            " μm is out of valid range [2, 1000] μm"
                        write(*,'(a)') "Valid wavelength range: 2.0 to 1000.0 μm"
                        if (present(is_success)) then
                            is_success = .false.
                            return
                        else
                            stop 1
                        end if
                    end if
                    if (this%wavelen_max(i) < 2.0_jprd .or. this%wavelen_max(i) > 1000.0_jprd) then
                        write(*,'(a,f8.2,a)') "ERROR: Wavelength max ", this%wavelen_max(i), &
                            " μm is out of valid range [2, 1000] μm"
                        write(*,'(a)') "Valid wavelength range: 2.0 to 1000.0 μm"
                        if (present(is_success)) then
                            is_success = .false.
                            return
                        else
                            stop 1
                        end if
                    end if
                end if
            end do
        end if

    end subroutine read_configuration

    !-----------------------------------------------------------------------
    ! Subroutine: initialize_hardware_config
    !
    ! Automatically detects hardware capabilities and optimizes configuration
    ! parameters for maximum performance on the current system.
    !-----------------------------------------------------------------------
    subroutine initialize_hardware_config(this)
        implicit none
        class(config_type), intent(inout) :: this
        
        ! Detect maximum available threads with SLURM/HPC awareness
        call detect_optimal_threads(this)
        
        ! Platform-specific optimizations
        call detect_platform_features(this)
        
        ! Auto-tune performance parameters based on hardware
        call optimize_for_hardware(this)
        
        ! Print detected configuration
        call this%print_config()
        
    end subroutine initialize_hardware_config
    
    !-----------------------------------------------------------------------
    ! Subroutine: detect_platform_features
    !
    ! Detects platform-specific features for optimization
    !-----------------------------------------------------------------------
    subroutine detect_platform_features(this)
        implicit none
        type(config_type), intent(inout) :: this
        character(len=256) :: sysctl_output, temp_file
        integer :: ios, unit_num
        logical :: file_exists
        
        ! Initialize with unknown
        this%cpu_architecture = "Unknown"
        unit_num = 99
        
#ifdef __APPLE__
        ! Apple platform detection
        this%cpu_architecture = "Apple"
        
        ! Try to detect Apple Silicon (M-series chips)
        open(unit=unit_num, file="/tmp/cpu_info.tmp", action="write", iostat=ios)
        if (ios == 0) then
            ! Use system command to get CPU info
            call execute_command_line("sysctl -n machdep.cpu.brand_string > /tmp/cpu_info.tmp 2>/dev/null", wait=.true.)
            close(unit_num)
            
            open(unit=unit_num, file="/tmp/cpu_info.tmp", action="read", iostat=ios)
            if (ios == 0) then
                read(unit_num, '(A)', iostat=ios) sysctl_output
                close(unit_num)
                if (ios == 0) then
                    if (index(sysctl_output, "Apple") > 0) then
                        this%is_apple_silicon = .true.
                        this%cpu_architecture = "Apple_Silicon"
                    end if
                end if
            end if
            
            ! Clean up temp file
            call execute_command_line("rm -f /tmp/cpu_info.tmp", wait=.false.)
        end if
#elif defined(__linux__) || defined(__linux) || defined(linux)
        ! Linux platform detection with multiple preprocessor checks
        this%cpu_architecture = "Linux"
#else
        ! Fallback: Try to detect Linux by checking for Linux-specific files
        inquire(file="/proc/version", exist=file_exists)
        if (file_exists) then
            this%cpu_architecture = "Linux"
        else
            inquire(file="/etc/os-release", exist=file_exists)
            if (file_exists) then
                this%cpu_architecture = "Linux"  
            else
                ! Try to detect using uname command
                temp_file = "/tmp/hpc_uname.tmp"
                call execute_command_line("uname -s > " // trim(temp_file) // " 2>/dev/null", wait=.true.)
                
                open(unit=unit_num, file=trim(temp_file), action="read", iostat=ios)
                if (ios == 0) then
                    read(unit_num, '(A)', iostat=ios) sysctl_output
                    close(unit_num)
                    if (ios == 0) then
                        if (index(sysctl_output, "Linux") > 0) then
                            this%cpu_architecture = "Linux"
                        endif
                    end if
                endif
                call execute_command_line("rm -f " // trim(temp_file), wait=.false.)
            endif
        endif
#endif
        
    end subroutine detect_platform_features
    
    !-----------------------------------------------------------------------
    ! Subroutine: detect_optimal_threads
    !
    ! Detects optimal thread count with SLURM/HPC environment awareness
    !-----------------------------------------------------------------------
    subroutine detect_optimal_threads(this)
        implicit none
        class(config_type), intent(inout) :: this
        character(len=128) :: env_value
        integer :: ios, detected_threads
        
        detected_threads = 0
        
        ! Priority 1: Check user-set OMP_NUM_THREADS
        call get_environment_variable("OMP_NUM_THREADS", env_value, status=ios)
        if (ios == 0 .and. len_trim(env_value) > 0) then
            read(env_value, *, iostat=ios) detected_threads
            if (ios == 0 .and. detected_threads > 0) then
                this%MAX_THREADS = detected_threads
                return
            endif
        endif
        
        ! Priority 2: Check SLURM environment variables
        call get_environment_variable("SLURM_CPUS_PER_TASK", env_value, status=ios)
        if (ios == 0 .and. len_trim(env_value) > 0) then
            read(env_value, *, iostat=ios) detected_threads
            if (ios == 0 .and. detected_threads > 0) then
                this%MAX_THREADS = detected_threads
                return
            endif
        endif
        
        call get_environment_variable("SLURM_NTASKS_PER_NODE", env_value, status=ios)
        if (ios == 0 .and. len_trim(env_value) > 0) then
            read(env_value, *, iostat=ios) detected_threads
            if (ios == 0 .and. detected_threads > 0) then
                this%MAX_THREADS = detected_threads
                return
            endif
        endif
        
        call get_environment_variable("SLURM_CPUS_ON_NODE", env_value, status=ios)
        if (ios == 0 .and. len_trim(env_value) > 0) then
            read(env_value, *, iostat=ios) detected_threads
            if (ios == 0 .and. detected_threads > 0) then
                this%MAX_THREADS = detected_threads
                return
            endif
        endif
        
        ! Priority 3: Check PBS/Torque environment variables
        call get_environment_variable("PBS_NUM_PPN", env_value, status=ios)
        if (ios == 0 .and. len_trim(env_value) > 0) then
            read(env_value, *, iostat=ios) detected_threads
            if (ios == 0 .and. detected_threads > 0) then
                this%MAX_THREADS = detected_threads
                return
            endif
        endif
        
        call get_environment_variable("NCPUS", env_value, status=ios)
        if (ios == 0 .and. len_trim(env_value) > 0) then
            read(env_value, *, iostat=ios) detected_threads
            if (ios == 0 .and. detected_threads > 0) then
                this%MAX_THREADS = detected_threads
                return
            endif
        endif
        
        ! Priority 4: Check LSF environment variables
        call get_environment_variable("LSB_DJOB_NUMPROC", env_value, status=ios)
        if (ios == 0 .and. len_trim(env_value) > 0) then
            read(env_value, *, iostat=ios) detected_threads
            if (ios == 0 .and. detected_threads > 0) then
                this%MAX_THREADS = detected_threads
                return
            endif
        endif
        
        ! Priority 5: Direct system detection for interactive HPC use
        call detect_system_cores_direct(detected_threads)
        if (detected_threads > 0) then
            this%MAX_THREADS = detected_threads
            return
        endif
        
        ! Priority 6: Fallback to OpenMP default
        this%MAX_THREADS = omp_get_max_threads()
        
        ! Ensure minimum of 1 thread
        if (this%MAX_THREADS < 1) then
            this%MAX_THREADS = 1
        endif
        
    end subroutine detect_optimal_threads
    
    !-----------------------------------------------------------------------
    ! Subroutine: detect_system_cores_direct
    !
    ! Directly detects CPU cores using system commands - critical for 
    ! interactive HPC usage where job scheduler variables are not set
    !-----------------------------------------------------------------------
    subroutine detect_system_cores_direct(detected_threads)
        implicit none
        integer, intent(out) :: detected_threads
        character(len=256) :: temp_file, temp_content
        character(len=128) :: env_value
        integer :: ios, unit_num, temp_threads
        
        detected_threads = 0
        unit_num = 98

#ifdef __APPLE__
        ! Apple macOS: Use sysctl for core count
        temp_file = "/tmp/macos_cores.tmp"
        call execute_command_line("sysctl -n hw.logicalcpu > " // trim(temp_file) // " 2>/dev/null", wait=.true.)
        
        open(unit=unit_num, file=trim(temp_file), action="read", iostat=ios)
        if (ios == 0) then
            read(unit_num, '(A)', iostat=ios) temp_content
            close(unit_num)
            if (ios == 0 .and. len_trim(temp_content) > 0) then
                read(temp_content, *, iostat=ios) temp_threads
                if (ios == 0 .and. temp_threads > 0) then
                    detected_threads = temp_threads
                    call execute_command_line("rm -f " // trim(temp_file), wait=.false.)
                    return
                endif
            endif
        endif
        call execute_command_line("rm -f " // trim(temp_file), wait=.false.)
#else
        ! Linux systems: Use nproc and /proc/cpuinfo
        
        ! Method 1: Try nproc command (most reliable on Linux)
        temp_file = "/tmp/hpc_nproc.tmp"
        call execute_command_line("nproc > " // trim(temp_file) // " 2>/dev/null", wait=.true.)
        
        open(unit=unit_num, file=trim(temp_file), action="read", iostat=ios)
        if (ios == 0) then
            read(unit_num, '(A)', iostat=ios) temp_content
            close(unit_num)
            if (ios == 0 .and. len_trim(temp_content) > 0) then
                read(temp_content, *, iostat=ios) temp_threads
                if (ios == 0 .and. temp_threads > 0) then
                    detected_threads = temp_threads
                    call execute_command_line("rm -f " // trim(temp_file), wait=.false.)
                    return
                endif
            endif
        endif
        call execute_command_line("rm -f " // trim(temp_file), wait=.false.)
        
        ! Method 2: Try parsing /proc/cpuinfo
        temp_file = "/tmp/hpc_cpuinfo.tmp"
        call execute_command_line("grep -c '^processor' /proc/cpuinfo > " // trim(temp_file) // " 2>/dev/null", wait=.true.)
        
        open(unit=unit_num, file=trim(temp_file), action="read", iostat=ios)
        if (ios == 0) then
            read(unit_num, '(A)', iostat=ios) temp_content
            close(unit_num)
            if (ios == 0 .and. len_trim(temp_content) > 0) then
                read(temp_content, *, iostat=ios) temp_threads
                if (ios == 0 .and. temp_threads > 0) then
                    detected_threads = temp_threads
                    call execute_command_line("rm -f " // trim(temp_file), wait=.false.)
                    return
                endif
            endif
        endif
        call execute_command_line("rm -f " // trim(temp_file), wait=.false.)
        
        ! Method 3: Try lscpu command
        temp_file = "/tmp/hpc_lscpu.tmp"
        call execute_command_line("lscpu -p | grep -v '^#' | wc -l > " // trim(temp_file) // " 2>/dev/null", wait=.true.)
        
        open(unit=unit_num, file=trim(temp_file), action="read", iostat=ios)
        if (ios == 0) then
            read(unit_num, '(A)', iostat=ios) temp_content
            close(unit_num)
            if (ios == 0 .and. len_trim(temp_content) > 0) then
                read(temp_content, *, iostat=ios) temp_threads
                if (ios == 0 .and. temp_threads > 0) then
                    detected_threads = temp_threads
                    call execute_command_line("rm -f " // trim(temp_file), wait=.false.)
                    return
                endif
            endif
        endif
        call execute_command_line("rm -f " // trim(temp_file), wait=.false.)
#endif
        
#ifndef __APPLE__
        ! Method 4: Fallback - try reading /proc/cpuinfo directly (Linux only)
        open(unit=unit_num, file="/proc/cpuinfo", action="read", iostat=ios)
        if (ios == 0) then
            temp_threads = 0
            do
                read(unit_num, '(A)', iostat=ios) temp_content
                if (ios /= 0) exit
                if (index(temp_content, "processor") == 1) then
                    temp_threads = temp_threads + 1
                endif
            end do
            close(unit_num)
            if (temp_threads > 0) then
                detected_threads = temp_threads
                return
            endif
        endif
#endif
        
        ! Method 5: Try environment variables that might be set on HPC systems
        call get_environment_variable("SLURM_JOB_CPUS_PER_NODE", env_value, status=ios)
        if (ios == 0 .and. len_trim(env_value) > 0) then
            read(env_value, *, iostat=ios) temp_threads
            if (ios == 0 .and. temp_threads > 0) then
                detected_threads = temp_threads
                return
            endif
        endif
        
        call get_environment_variable("NSLOTS", env_value, status=ios)
        if (ios == 0 .and. len_trim(env_value) > 0) then
            read(env_value, *, iostat=ios) temp_threads
            if (ios == 0 .and. temp_threads > 0) then
                detected_threads = temp_threads
                return
            endif
        endif
        
    end subroutine detect_system_cores_direct
    
    !-----------------------------------------------------------------------
    ! Subroutine: optimize_for_hardware
    !
    ! Optimizes configuration parameters based on detected hardware
    !-----------------------------------------------------------------------
    subroutine optimize_for_hardware(this)
        implicit none
        class(config_type), intent(inout) :: this
        integer(kind=jpim) :: problem_size, work_per_thread
        
        if (this%is_apple_silicon) then
            ! Apple M-series optimization (enhanced for performance)
            this%L1_CACHE_SIZE        = 262144    ! 256KB L1 cache per core
            this%CACHE_LINE_SIZE      = 128       ! M3 128-byte cache line
            this%VECTOR_LENGTH        = 32        ! M3 SIMD width (256-bit)
            this%L2_BLOCK_SIZE        = 2048      ! Reduced for better cache fit
            this%MIN_CHUNK_SIZE       = 64        ! Optimized for vectorization
            this%MAX_CHUNK_SIZE       = 1024      ! Reduced to prevent cache thrashing
            this%MIN_PARALLEL_SIZE    = 256       ! Lower threshold for better utilization
            this%MIN_WORK_PER_THREAD  = 128       ! Balanced work distribution
        else
            ! Enhanced optimization for Intel/AMD processors
            this%L1_CACHE_SIZE        = 32768     ! 32KB typical L1
            this%CACHE_LINE_SIZE      = 64        ! Standard cache line
            this%VECTOR_LENGTH        = 8         ! AVX2 width
            this%L2_BLOCK_SIZE        = 512       ! Optimized for cache efficiency
            this%MIN_CHUNK_SIZE       = 32        ! Smaller chunks for better balance
            this%MAX_CHUNK_SIZE       = 512       ! Reduced to improve cache utilization
            this%MIN_PARALLEL_SIZE    = 100       ! Lower threshold for better utilization
            this%MIN_WORK_PER_THREAD  = 64        ! Optimized work distribution
        end if
        
        ! Scale parameters based on thread count for optimal performance
        if (this%MAX_THREADS >= 16) then
            this%MIN_PARALLEL_SIZE = this%MIN_PARALLEL_SIZE / 2  ! More aggressive parallelization
            this%MAX_CHUNK_SIZE = this%MAX_CHUNK_SIZE * 2
        elseif (this%MAX_THREADS >= 8) then
            this%MIN_PARALLEL_SIZE = this%MIN_PARALLEL_SIZE
            this%MAX_CHUNK_SIZE = int(real(this%MAX_CHUNK_SIZE) * 1.5)
        elseif (this%MAX_THREADS >= 4) then
            this%MIN_PARALLEL_SIZE = this%MIN_PARALLEL_SIZE * 2
        end if
        
        ! NEW: Optimize based on problem size (n_angles × n_winds)
        ! This ensures performance parameters adapt to actual workload
        if (this%n_angles > 0 .and. this%n_winds > 0) then
            problem_size = this%n_angles * this%n_winds
            work_per_thread = problem_size / max(1, this%MAX_THREADS)
            
            ! Adjust chunk sizes based on problem dimensions
            if (problem_size > 10000) then
                ! Large problem: increase chunk sizes for better cache utilization
                this%MAX_CHUNK_SIZE = min(2048, this%MAX_CHUNK_SIZE * 2)
                this%MIN_CHUNK_SIZE = min(256, this%MIN_CHUNK_SIZE * 2)
            elseif (problem_size < 100) then
                ! Small problem: reduce parallel overhead
                this%MIN_PARALLEL_SIZE = problem_size
                this%MIN_CHUNK_SIZE = max(1, problem_size / this%MAX_THREADS)
                this%MAX_CHUNK_SIZE = problem_size
            end if
            
            ! Adjust work distribution based on actual work per thread
            if (work_per_thread < 10) then
                ! Very little work per thread: disable parallelization
                this%MIN_PARALLEL_SIZE = problem_size + 1  ! Force serial execution
            elseif (work_per_thread < 50) then
                ! Limited work: reduce chunk sizes for better load balancing
                this%MIN_CHUNK_SIZE = max(1, work_per_thread / 4)
                this%MAX_CHUNK_SIZE = work_per_thread
            end if
            
            ! Special handling for extreme array dimensions
            if (this%n_angles > 100 .or. this%n_winds > 100) then
                ! Large arrays: optimize for memory bandwidth
                this%use_cache_blocking = .true.
                this%L2_BLOCK_SIZE = min(this%L2_BLOCK_SIZE, 256)
            end if
        end if
        
        ! Enable optimizations based on hardware
        this%use_guided_scheduling = .true.  ! Always use guided scheduling
        this%use_simd_optimization = .true.  ! Always enable SIMD
        this%use_cache_blocking = this%use_cache_blocking .or. (this%MAX_THREADS >= 4)  ! Enable for multi-core
        
    end subroutine optimize_for_hardware
    
    !-----------------------------------------------------------------------
    ! Subroutine: print_configuration
    !
    ! Prints the current configuration for debugging and verification
    !-----------------------------------------------------------------------
    subroutine print_configuration(this)
        implicit none
        class(config_type), intent(in) :: this
        
        write(*,'(A)') "======================================"
        write(*,'(A)') "   Hardware Configuration Detected   "
        write(*,'(A)') "======================================"
        write(*,'(A,I0)') "Maximum Threads:      ", this%MAX_THREADS
        write(*,'(A,A)')  "CPU Architecture:     ", trim(this%cpu_architecture)
        write(*,'(A,L1)') "Apple Silicon:        ", this%is_apple_silicon
        write(*,'(A)')    "--------------------------------------"
        write(*,'(A,I0)') "L1 Cache Size:        ", this%L1_CACHE_SIZE
        write(*,'(A,I0)') "Cache Line Size:      ", this%CACHE_LINE_SIZE
        write(*,'(A,I0)') "Vector Length:        ", this%VECTOR_LENGTH
        write(*,'(A,I0)') "L2 Block Size:        ", this%L2_BLOCK_SIZE
        write(*,'(A)')    "--------------------------------------"
        write(*,'(A,I0)') "Min Chunk Size:       ", this%MIN_CHUNK_SIZE
        write(*,'(A,I0)') "Max Chunk Size:       ", this%MAX_CHUNK_SIZE
        write(*,'(A,I0)') "Min Parallel Size:    ", this%MIN_PARALLEL_SIZE
        write(*,'(A,I0)') "Min Work Per Thread:  ", this%MIN_WORK_PER_THREAD
        write(*,'(A)')    "--------------------------------------"
        write(*,'(A,L1)') "Guided Scheduling:    ", this%use_guided_scheduling
        write(*,'(A,L1)') "SIMD Optimization:    ", this%use_simd_optimization
        write(*,'(A,L1)') "Cache Blocking:       ", this%use_cache_blocking
        write(*,'(A)') "======================================"
        write(*,'(A)') ""
        
    end subroutine print_configuration

end module configuration
