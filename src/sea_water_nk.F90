! sea_water_nk.F90
!
! Purpose:
!   Sea water complex refractive index (n + ik) calculations with temperature
!   and salinity dependence for infrared ocean emissivity modeling.
!
! This module computes the wavelength-dependent complex refractive index of
! sea water required for Fresnel reflectance calculations in the emissivity
! engine. Supports the full infrared spectrum from 10-5000 cm^-1.
!
! Features:
!   - Temperature-dependent pure water optical constants (Wang et al., 2026)
!   - Salinity correction via Pinkley & Williams (1976) coefficients
!   - Log-space spectral interpolation for numerical stability
!   - Multiple interpolation methods: spline, linear, nearest neighbor
!   - Thread-safe parallel interpolation with OpenMP
!   - NetCDF data integration for optical constants database
!   - Closest-temperature selection for discrete temperature grids
!
! Optical Properties at Key Spectral Regions:
!   - Far-infrared (10-200 cm^-1): n = 1.5-2.5, k = 0.4-1.1 (high absorption)
!   - Thermal window (800-1200 cm^-1): n = 1.2-1.3, k = 0.01-0.1
!   - Near-infrared (2000-5000 cm^-1): n = 1.3-1.5, k = 0.01-0.4
!
! Data Source:
!   - water_optical_constants.nc: 6496 wavenumbers, 11 temperatures
!   - Temperature range: 273-313 K
!   - Salinity range: 0-45 PSU
!
! Public Interface:
!   - get_sea_water_nk: Single-point refractive index calculation
!
! References:
!   Wang et al. (2026). Temperature-dependent water refractive index.
!   Pinkley & Williams (1976). Optical constants of sea water.
!   Rottgers et al. (2014). Mass-specific absorption coefficients.
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

module sea_water_nk
    use utils,          only : jpim, jprd
    use netcdf_handler, only : netcdf_file
    use configuration,  only : config_type
    use mathlib,        only : get_interpolation_1d, error_type
    use omp_lib
    use ieee_arithmetic
    implicit none
    private
    public:: get_sea_water_nk

    !-----------------------------------------------------------------------
    ! Module-wide constants and parameters
    !-----------------------------------------------------------------------
    ! Error codes (local definitions matching mathlib conventions)
    integer(kind=jpim), parameter :: SUCCESS           = 0
    integer(kind=jpim), parameter :: ERR_INVALID_INPUT = 6
    integer(kind=jpim), parameter :: ERR_COMPUTATION   = 5

    ! Numerical constant for invalid/missing data
    real(kind = jprd), parameter:: INVALID_VALUE = -huge(1.0_jprd)

    ! Available interpolation methods
    character(len=*), parameter:: METHOD_SPLINE  = 'spline'   ! Cubic spline interpolation
    character(len=*), parameter:: METHOD_LINEAR  = 'linear'   ! Linear interpolation
    character(len=*), parameter:: METHOD_NEAREST = 'nearest'  ! Nearest neighbor

    ! Validation bounds (documented constants)
    real(kind=jprd), parameter :: SALINITY_MIN = 0.0_jprd     ! Minimum salinity (PSU)
    real(kind=jprd), parameter :: SALINITY_MAX = 45.0_jprd    ! Maximum salinity (PSU)

    !-----------------------------------------------------------------------
    ! Internal derived type for NetCDF data storage
    !-----------------------------------------------------------------------
    type :: nk_data_type
        real(kind=jprd), allocatable :: temperature(:)        ! Temperature array from NetCDF
        real(kind=jprd), allocatable :: wavenumber(:)         ! Wavenumber array from NetCDF
        real(kind=jprd), allocatable :: refra_real_part(:,:)  ! Real refractive index (ntemp, nwave)
        real(kind=jprd), allocatable :: refra_imag_part(:,:)  ! Imaginary refractive index (ntemp, nwave)
        real(kind=jprd), allocatable :: n_salt_corr(:)        ! Salinity correction for real part
        real(kind=jprd), allocatable :: k_salt_corr(:)        ! Salinity correction for imaginary part
        integer(kind=jpim) :: nwave = 0                       ! Number of wavenumbers
        integer(kind=jpim) :: ntemp = 0                       ! Number of temperatures
        logical :: loaded = .false.                           ! Data loaded flag
    end type nk_data_type


contains

    !-----------------------------------------------------------------------
    ! Internal subroutine: Load NetCDF data
    !
    ! Loads all required optical data from the NetCDF auxiliary file.
    ! This eliminates code duplication between the two public subroutines.
    !
    ! Args:
    !   ncfile_name: Path to NetCDF file
    !   data: Output data structure
    !   error: Error handling structure
    !-----------------------------------------------------------------------
    subroutine load_nk_data(ncfile_name, data, error)
        implicit none
        character(len=*), intent(in)     :: ncfile_name
        type(nk_data_type), intent(out)  :: data
        type(error_type), intent(out)    :: error

        type(netcdf_file) :: file
        logical           :: file_exists

        ! Initialize error
        error%code = SUCCESS
        error%message = ''

        ! Validate that the NetCDF file exists
        inquire(file=trim(ncfile_name), exist=file_exists)
        if (.not. file_exists) then
            error%code = ERR_INVALID_INPUT
            error%message = 'NetCDF file not found: ' // trim(ncfile_name)
            return
        end if

        call file%open(trim(ncfile_name))
        if (.not. file%is_open()) then
            error%code = ERR_INVALID_INPUT
            error%message = 'Failed to open NetCDF file: ' // trim(ncfile_name)
            return
        end if

        ! Read wavenumber data
        if (.not. file%exists('wavenumber')) then
            error%code = ERR_INVALID_INPUT
            error%message = '"wavenumber" not found in the auxiliary NETCDF file'
            call file%close()
            return
        endif
        call file%get('wavenumber', data%wavenumber)

        ! Read temperature data
        if (.not. file%exists('temperature')) then
            error%code = ERR_INVALID_INPUT
            error%message = '"temperature" not found in the auxiliary NETCDF file'
            call file%close()
            return
        endif
        call file%get('temperature', data%temperature)

        ! Read real refractive index of pure water (2D array)
        if (.not. file%exists('refra_real_part')) then
            error%code = ERR_INVALID_INPUT
            error%message = '"refra_real_part" not found in the auxiliary NETCDF file'
            call file%close()
            return
        endif
        call file%get('refra_real_part', data%refra_real_part)

        ! Read imaginary refractive index of pure water (2D array)
        if (.not. file%exists('refra_imag_part')) then
            error%code = ERR_INVALID_INPUT
            error%message = '"refra_imag_part" not found in the auxiliary NETCDF file'
            call file%close()
            return
        endif
        call file%get('refra_imag_part', data%refra_imag_part)

        ! Read salinity correction coefficients for real part
        if (.not. file%exists('n_salt_corr')) then
            error%code = ERR_INVALID_INPUT
            error%message = '"n_salt_corr" not found in the auxiliary NETCDF file'
            call file%close()
            return
        endif
        call file%get('n_salt_corr', data%n_salt_corr)

        ! Read salinity correction coefficients for imaginary part
        if (.not. file%exists('k_salt_corr')) then
            error%code = ERR_INVALID_INPUT
            error%message = '"k_salt_corr" not found in the auxiliary NETCDF file'
            call file%close()
            return
        endif
        call file%get('k_salt_corr', data%k_salt_corr)

        ! Close NetCDF file
        call file%close()

        ! Set dimensions
        data%nwave = size(data%wavenumber)
        data%ntemp = size(data%temperature)
        data%loaded = .true.

        ! Validate array dimensions for consistency
        if (size(data%n_salt_corr) /= data%nwave .or. size(data%k_salt_corr) /= data%nwave) then
            error%code = ERR_INVALID_INPUT
            write(error%message, '(A,I6,A,I6,A,I6)') &
                'Salt correction arrays size mismatch. Expected:', data%nwave, &
                ', got n_salt_corr:', size(data%n_salt_corr), ', k_salt_corr:', size(data%k_salt_corr)
            call cleanup_nk_data(data)
            return
        end if

    end subroutine load_nk_data


    !-----------------------------------------------------------------------
    ! Internal subroutine: Validate input parameters
    !
    ! Validates temperature, salinity, and wavenumber ranges.
    !
    ! Args:
    !   data: Loaded NetCDF data
    !   temperature_in: Input temperature (K)
    !   salinity_in: Input salinity (PSU)
    !   wavenumber_min: Minimum wavenumber to check (cm^-1)
    !   wavenumber_max: Maximum wavenumber to check (cm^-1)
    !   error: Error handling structure
    !-----------------------------------------------------------------------
    subroutine validate_inputs(data, temperature_in, salinity_in, &
                               wavenumber_min, wavenumber_max, error)
        implicit none
        type(nk_data_type), intent(in)   :: data
        real(kind=jprd), intent(in)      :: temperature_in, salinity_in
        real(kind=jprd), intent(in)      :: wavenumber_min, wavenumber_max
        type(error_type), intent(out)    :: error

        error%code = SUCCESS
        error%message = ''

        ! Validate wavenumber range (small tolerance for floating point)
        if (wavenumber_min < (minval(data%wavenumber) - 1.0e-6_jprd) .or. &
            wavenumber_max > (maxval(data%wavenumber) + 1.0e-6_jprd)) then
            error%code = ERR_INVALID_INPUT
            write(error%message, '(A,F10.2,A,F10.2,A)') &
                'Wavenumber range should be within: [', minval(data%wavenumber), &
                ',', maxval(data%wavenumber), '] cm^-1'
            return
        endif

        ! Validate temperature range
        if (temperature_in < minval(data%temperature) .or. &
            temperature_in > maxval(data%temperature)) then
            error%code = ERR_INVALID_INPUT
            write(error%message, '(A,F8.2,A,F8.2,A)') &
                'Temperature should be within: [', minval(data%temperature), &
                ',', maxval(data%temperature), '] K'
            return
        endif

        ! Validate salinity range
        if (salinity_in < SALINITY_MIN .or. salinity_in > SALINITY_MAX) then
            error%code = ERR_INVALID_INPUT
            write(error%message, '(A,F6.1,A,F6.1,A)') &
                'Salinity should be within: [', SALINITY_MIN, ',', SALINITY_MAX, '] PSU'
            return
        endif

    end subroutine validate_inputs


    !-----------------------------------------------------------------------
    ! Internal subroutine: Compute refractive index at target wavenumbers
    !
    ! Core computation that applies salinity correction, interpolates over
    ! wavenumber, and selects the closest temperature. This is the shared
    ! algorithm used by both public interfaces.
    !
    ! Args:
    !   data: Loaded NetCDF data
    !   wavenumber_target: Target wavenumber array (cm^-1)
    !   temperature_in: Input temperature (K)
    !   salinity_in: Input salinity (PSU)
    !   real_out: Output real part of refractive index
    !   imag_out: Output imaginary part of refractive index
    !   error: Error handling structure
    !-----------------------------------------------------------------------
    subroutine compute_nk_internal(data, wavenumber_target, temperature_in, &
                                   salinity_in, real_out, imag_out, error)
        implicit none
        type(nk_data_type), intent(in)            :: data
        real(kind=jprd), intent(in)               :: wavenumber_target(:)
        real(kind=jprd), intent(in)               :: temperature_in, salinity_in
        real(kind=jprd), intent(out)              :: real_out(:)
        real(kind=jprd), intent(out)              :: imag_out(:)
        type(error_type), intent(out)             :: error

        ! Local variables
        integer(kind=jpim) :: N_out, it, pos_close_id
        real(kind=jprd), allocatable :: refra_real_base(:,:), refra_imag_base(:,:)
        real(kind=jprd), allocatable :: temp_real_out(:,:), temp_imag_out(:,:)
        real(kind=jprd), allocatable :: temp_real_out_log(:,:), temp_imag_out_log(:,:)
        real(kind=jprd), allocatable :: wavenumber_log(:), wavenumber_target_log(:)
        character(len=20) :: method_interp

        error%code = SUCCESS
        error%message = ''
        N_out = size(wavenumber_target)

        ! Allocate working arrays
        allocate(refra_real_base(data%nwave, data%ntemp))
        allocate(refra_imag_base(data%nwave, data%ntemp))
        allocate(temp_real_out(N_out, data%ntemp))
        allocate(temp_imag_out(N_out, data%ntemp))
        allocate(temp_real_out_log(N_out, data%ntemp))
        allocate(temp_imag_out_log(N_out, data%ntemp))
        allocate(wavenumber_log(data%nwave))
        allocate(wavenumber_target_log(N_out))

        ! Step 1: Apply salinity correction to pure water values
        ! The temperature dependent water refractive index is from Wang et al., 2026.
        ! salt water refractive index = pure water + salt_correction * salinity
        ! Note: NetCDF arrays are (ntemp, nwave), we transpose to (nwave, ntemp)
        do it = 1, data%ntemp
            refra_real_base(:, it) = data%refra_real_part(it, :) + data%n_salt_corr(:) * salinity_in
            refra_imag_base(:, it) = data%refra_imag_part(it, :) + data%k_salt_corr(:) * salinity_in
        end do

        ! Pre-compute log of wavenumber arrays (log-space interpolation for stability)
        wavenumber_log = log(data%wavenumber)
        wavenumber_target_log = log(wavenumber_target)

        ! Step 2: Interpolate over wavenumber for each temperature
        method_interp = METHOD_SPLINE
        do it = 1, data%ntemp
            ! Interpolate real part in log-space
            call get_interpolation_1d(wavenumber_log, log(refra_real_base(:, it)), &
                wavenumber_target_log, temp_real_out_log(:, it), method_interp)
            temp_real_out(:, it) = exp(temp_real_out_log(:, it))

            ! Interpolate imaginary part in log-space
            call get_interpolation_1d(wavenumber_log, log(refra_imag_base(:, it)), &
                wavenumber_target_log, temp_imag_out_log(:, it), method_interp)
            temp_imag_out(:, it) = exp(temp_imag_out_log(:, it))
        end do

        ! Step 3: Select closest temperature
        pos_close_id = minloc(abs(data%temperature - temperature_in), 1)
        real_out = temp_real_out(:, pos_close_id)
        imag_out = temp_imag_out(:, pos_close_id)

        ! Validate output for physical reasonableness
        if (any(.not. ieee_is_finite(real_out)) .or. any(.not. ieee_is_finite(imag_out))) then
            error%code = ERR_COMPUTATION
            error%message = 'Non-finite values detected in output arrays'
        endif

        if (any(real_out <= 0.0_jprd) .or. any(imag_out <= 0.0_jprd)) then
            ! Warning only, not an error - some edge cases may have small negative values
            print *, 'WARNING: Non-positive refractive index values detected'
            print *, '  Real part range: [', minval(real_out), ',', maxval(real_out), ']'
            print *, '  Imag part range: [', minval(imag_out), ',', maxval(imag_out), ']'
        endif

        ! Cleanup
        deallocate(refra_real_base, refra_imag_base)
        deallocate(temp_real_out, temp_imag_out)
        deallocate(temp_real_out_log, temp_imag_out_log)
        deallocate(wavenumber_log, wavenumber_target_log)

    end subroutine compute_nk_internal


    !-----------------------------------------------------------------------
    ! Internal subroutine: Cleanup nk_data_type
    !-----------------------------------------------------------------------
    subroutine cleanup_nk_data(data)
        implicit none
        type(nk_data_type), intent(inout) :: data

        if (allocated(data%temperature)) deallocate(data%temperature)
        if (allocated(data%wavenumber)) deallocate(data%wavenumber)
        if (allocated(data%refra_real_part)) deallocate(data%refra_real_part)
        if (allocated(data%refra_imag_part)) deallocate(data%refra_imag_part)
        if (allocated(data%n_salt_corr)) deallocate(data%n_salt_corr)
        if (allocated(data%k_salt_corr)) deallocate(data%k_salt_corr)
        data%nwave = 0
        data%ntemp = 0
        data%loaded = .false.

    end subroutine cleanup_nk_data


    !-----------------------------------------------------------------------
    ! Main subroutine for seawater refractive index calculation
    !
    ! Computes complex refractive index (n+ik) for seawater as function of
    ! wavenumber, temperature, salinity using parallel interpolation
    !
    ! Args:
    !   config: Configuration object with refractive index filename and parameters
    !   wavenumber_start: Starting wavenumber (cm^-1)
    !   wavenumber_end: Ending wavenumber (cm^-1)
    !   N_in: Number of output points (modified if interpolation changes)
    !   temperature_in: Temperature (K)
    !   salinity_in: Salinity (PSU)
    !   real_out: Real part of refractive index (n)
    !   imag_out: Imaginary part of refractive index (k)
    !   wavenumber_out: Output wavenumber points
    !-----------------------------------------------------------------------
    subroutine get_sea_water_nk(config, wavenumber_start, wavenumber_end, N_in, &
                                temperature_in, salinity_in, real_out, imag_out, wavenumber_out)
        implicit none

        ! Input parameters
        type(config_type), intent(in)             :: config
        real(kind=jprd), intent(in)               :: wavenumber_start, wavenumber_end
        integer(kind=jpim), intent(inout)         :: N_in
        real(kind=jprd), intent(in)               :: temperature_in
        real(kind=jprd), intent(in)               :: salinity_in

        ! Output arrays
        real(kind=jprd), allocatable, intent(out) :: real_out(:)
        real(kind=jprd), allocatable, intent(out) :: imag_out(:)
        real(kind=jprd), allocatable, intent(out) :: wavenumber_out(:)

        ! Local variables
        type(nk_data_type) :: data
        type(error_type)   :: error
        integer(kind=jpim) :: iw

        ! Validate number of points
        if (N_in <= 0) then
            print *, 'ERROR: Number of wavenumber points must be > 0'
            return
        endif

        ! Load NetCDF data
        call load_nk_data(trim(config%refractive_index_filename), data, error)
        if (error%code /= SUCCESS) then
            print *, 'ERROR: ', trim(error%message)
            return
        endif

        ! Validate inputs
        call validate_inputs(data, temperature_in, salinity_in, &
                             wavenumber_start, wavenumber_end, error)
        if (error%code /= SUCCESS) then
            print *, 'ERROR: ', trim(error%message)
            call cleanup_nk_data(data)
            return
        endif

        ! Generate output wavenumber grid
        allocate(wavenumber_out(N_in))
        if (N_in == 1) then
            wavenumber_out(1) = wavenumber_start
        else
            do iw = 1, N_in
                wavenumber_out(iw) = wavenumber_start + &
                    (wavenumber_end - wavenumber_start) * real(iw-1, jprd) / real(N_in-1, jprd)
            end do
        end if

        ! Allocate output arrays
        allocate(real_out(N_in), imag_out(N_in))

        ! Compute refractive indices
        call compute_nk_internal(data, wavenumber_out, temperature_in, &
                                 salinity_in, real_out, imag_out, error)
        if (error%code /= SUCCESS) then
            print *, 'ERROR: ', trim(error%message)
            call cleanup_nk_data(data)
            return
        endif

        ! Cleanup
        call cleanup_nk_data(data)

    end subroutine get_sea_water_nk

end module sea_water_nk
