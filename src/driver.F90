! driver.F90 - Main Driver for IOSEE Ocean Emissivity Package
!
! Purpose:
!   This is the main program for running the
!   Infrared Ocean Surface Effective Emissivity (IOSEE) package with support for:
!   - Spectral mode: High-resolution wavenumber calculations (10-5000 cm^-1)
!   - Narrow-band mode: Band-averaged emissivity with Planck integration
!   - Polarized mode: Separate V and H polarization emissivity output
!
!   Reads configuration and input data, performs emissivity calculations using
!   Wu & Smith (1997) and Masuda (2006) physics, outputs results to NetCDF.
!
! Physical Models:
!   - Wu & Smith (1997): Rough sea surface emissivity via Cox-Munk integration
!   - Masuda (2006): Multiple reflection effects between wave facets
!   - Nalli et al. (2008, 2023): Effective view angle corrections
!     SST-dependent 3D LUT (71 angles x 19 winds x 21 SSTs), PCHIP-interpolated
!     along SST at runtime to produce a 2D (angle, wind) slice
!   - Henderson et al. (2003): Polarized emissivity modeling
!
! Performance Optimizations:
!   - Single NetCDF read for optical constants (eliminates redundant I/O)
!   - Pre-allocated arrays outside processing loops
!   - OpenMP parallelization with COLLAPSE(3) and guided scheduling
!   - SIMD vectorization with 32-byte alignment
!   - Cross-platform: Apple Silicon M-series and Linux HPC systems
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
! Command Line Arguments (flexible key=value format, order independent):
!   config_file              : Configuration file (namelist format, required)
!   output=<file>            : Output file path (NetCDF, optional override)
!
! Configuration Options:
!   mode = 'spectral'        : Continuous wavenumber calculations
!   mode = 'narrow-band'     : Band-averaged with Planck function weighting
!   polarization = 'enable'  : Output V and H polarization separately
!

program driver

    use utils,              only: jpim, jprd
    use netcdf_handler,     only: netcdf_file
    use configuration,      only: config_type
    use mathlib,            only: initialize_parallel_env, error_type, effective_view_angle, &
                                  smooth_effective_angle_lut, interpolate_lut_sst
    use ocean_emissivity,   only: get_emissivity_optimized, &
                                  get_polarized_emissivity_optimized, &
                                  initialize_workspace_pool_emiss
    use sea_water_nk,       only: get_sea_water_nk
    use save_output,        only: save_spectral_emissivity, save_narrow_band_emissivity, save_spectral_emissivity_wl, &
                                  save_spectral_polarized_emissivity, save_spectral_polarized_emissivity_wl
    use omp_lib
    implicit none

    ! Local declarations
    type(config_type)   :: config        ! Configuration object
    character(len = 512)  :: config_file   ! Configuration file path
    character(len = 512)  :: output_file   ! NetCDF output file path
    character(len = 512)  :: arg_string    ! Command line argument
    character(len = 512)  :: arg_lower     ! Lowercase version of argument
    integer(kind = jpim)  :: iarg          ! Command line argument loop index
    logical            :: success        ! Configuration read success flag
    integer(kind = jpim):: istatus        ! Command line argument status
    integer(kind = jpim):: num_threads    ! Number of OpenMP threads
    integer(kind = jpim):: ierr           ! General error status
    character(len = 8)   :: thread_str     ! Thread number string conversion
    integer(kind = jpim):: n_angles       ! Number of viewing angles
    integer(kind = jpim):: n_winds        ! Number of wind speeds
    integer(kind = jpim):: n_wavenum      ! Number of wavenumbers
    integer(kind = jpim):: i, j, k        ! Loop indices
    real(kind = jprd)    :: temp_in        ! Input temperature
    real(kind = jprd)    :: salinity_in    ! Input salinity
    real(kind = jprd)    :: wavenum_start  ! Starting wavenumber
    real(kind = jprd)    :: wavenum_end    ! Ending wavenumber
    type(error_type)   :: err_oe         ! Ocean emissivity error
    logical            :: is_valid       ! Input validation flag
    character(len = 256):: ncfile_name   ! Netcdf file

    ! Arrays for spectral calculations
    real(kind = jprd), allocatable:: wavenum(:)    ! Wavenumber array
    real(kind = jprd), allocatable:: wavelen(:)    ! Wavelength array
    real(kind = jprd), allocatable:: real_n(:)     ! Real refractive index
    real(kind = jprd), allocatable:: imag_n(:)     ! Imaginary refractive index
    complex(kind = jprd), allocatable:: refm(:)    ! Complex refractive index

    ! Arrays for angles and winds
    real(kind = jprd), allocatable:: angles(:)     ! Viewing angles
    real(kind = jprd), allocatable:: winds(:)      ! Wind speeds
    real(kind = jprd), allocatable:: effective_angles(:,:)  ! Effective viewing angles (n_angles, n_winds)

    ! Arrays for results
    real(kind = jprd), allocatable:: emiss(:,:,:)        ! Final emissivity
    real(kind = jprd), allocatable:: emiss_noref(:,:,:)  ! Emissivity without reflection
    real(kind = jprd), allocatable:: refl_cb(:,:,:)      ! Cox-Munk reflectance

    ! For timing and performance
    integer(kind = jpim):: start_clock, end_clock, clock_rate
    real(kind = jprd)    :: elapsed_time
    real(kind = jprd)    :: performance_rate

    ! For effective view angle calculation from NetCDF
    type(netcdf_file)  :: file           ! NetCDF file handle
    logical            :: file_exists    ! File existence check
    real(kind = jprd), allocatable:: reference_angle(:)                     ! Reference angles from NetCDF
    real(kind = jprd), allocatable:: reference_wind_speed(:)               ! Reference wind speeds from NetCDF
    real(kind = jprd), allocatable:: effective_view_angle_lut(:,:)       ! Lookup table
    real(kind = jprd), allocatable:: effective_view_angle_lut_3d(:,:,:) ! 3D LUT from file
    real(kind = jprd), allocatable:: reference_sst(:)                   ! SST coordinate from LUT

    ! For normal view angle and angle blending at 70° boundary
    real(kind = jprd), allocatable:: reference_angle_nor(:)               ! Reference angles from normal LUT
    real(kind = jprd), allocatable:: reference_wind_speed_nor(:)         ! Reference wind speeds from normal LUT
    real(kind = jprd), allocatable:: normal_view_angle_lut(:,:)          ! Normal view angle lookup table
    real(kind = jprd), allocatable:: angles_eff_calc(:)                  ! Angles for effective calc (<=70° + 70° bridge)
    real(kind = jprd), allocatable:: angles_nor_calc(:)                  ! Angles for normal calc (70° + >70°)
    real(kind = jprd), allocatable:: eff_angles_eff(:,:)                 ! Effective angles from eff LUT
    real(kind = jprd), allocatable:: eff_angles_nor(:,:)                 ! Effective angles from normal LUT
    logical            :: has_angles_gt70   ! Whether any user angles > 70°
    logical            :: has_70_exact      ! Whether 70 deg is exactly in the user angle list
    integer(kind = jpim):: n_angles_le70    ! Count of user angles <= 70°
    integer(kind = jpim):: n_angles_gt70    ! Count of user angles > 70°
    integer(kind = jpim):: n_eff_calc       ! Size of angles_eff_calc
    integer(kind = jpim):: n_nor_calc       ! Size of angles_nor_calc
    integer(kind = jpim):: idx_70_eff       ! Index of 70° in angles_eff_calc

    !----------------------------------------------------------------------
    ! 0) Parse command line arguments and read configuration
    !----------------------------------------------------------------------
    ! Parse command line: ./iosee config.file [output=<output.nc>]
    if (command_argument_count() < 1) then
        write(*,'(a)') "Usage: ./iosee config.file [output=<filename.nc>]"
        write(*,'(a)') "  config.file        - Configuration file path (required)"
        write(*,'(a)') "  output=<file.nc>   - Output NetCDF file (optional, default: emissivity_output.nc)"
        stop 1
    end if

    ! Get config file (first argument)
    call get_command_argument(1, config_file, status = istatus)
    if (istatus /= 0) then
        write(*,'(a)') "Error: Could not read config file argument"
        stop 1
    end if

    ! Get output file (search all arguments for output= prefix, case-insensitive)
    ! Initialize output_file to empty to detect command-line override
    output_file = ""

    ! Search all arguments for output= prefix (case-insensitive)
    do iarg = 2, command_argument_count()
        call get_command_argument(iarg, arg_string, status = istatus)
        if (istatus == 0) then
            ! Convert to lowercase for comparison
            arg_lower = arg_string
            call to_lowercase(arg_lower)
            arg_lower = adjustl(arg_lower)
            if (index(arg_lower, "output=") == 1) then
                ! Extract filename preserving original case, skip "output="
                arg_string = adjustl(arg_string)
                output_file = trim(adjustl(arg_string(8:)))
                exit
            end if
        end if
    end do

    write(*,'(a)') "================================================="
    write(*,'(a)') "   Infrared Ocean Surface Effective Emissivity (IOSEE) "
    write(*,'(a)') "================================================="
    write(*,'(a, a)') "Config file:  ", trim(config_file)

    ! Initialize hardware configuration and auto-detection
    call config%initialize()

    ! Read configuration (which may override auto-detected values)
    call config%read(config_file, is_success = success)
    if (.not. success) then
        write(*,'(a, a)') "Error: Failed to read config file: ", trim(config_file)
        stop 1
    end if

    ! Apply output file precedence: command-line > config file > default
    if (len_trim(output_file) == 0) then
        if (len_trim(config%output_file) > 0) then
            output_file = trim(config%output_file)
        else
            output_file = "emissivity_output.nc"
        end if
    end if

    write(*,'(a, a)') "Output file:  ", trim(output_file)
    write(*,'(a)') ""

    ! Re-optimize performance parameters based on actual problem size from config
    ! This second optimization pass considers n_angles and n_winds
    call config%optimize_for_hardware()

    !----------------------------------------------------------------------
    ! 1) Initialize parallel environment and workspace pools
    !----------------------------------------------------------------------
    ! Set thread count from detected hardware configuration
    num_threads = config%MAX_THREADS
    call initialize_parallel_env(num_threads)
    print *, "Initialized parallel environment with", num_threads, "OpenMP threads."

    call initialize_workspace_pool_emiss(err_oe)
    if (err_oe%code /= 0) then
        print *, "Error initializing ocean emissivity pool:", err_oe%message
        stop
    end if
    print *, "Ocean emissivity workspace pool initialized."

    !----------------------------------------------------------------------
    ! 2) Set up calculation parameters (now from config file)
    !----------------------------------------------------------------------
    ! Read parameters from config
    n_angles = config%n_angles
    n_winds  = config%n_winds
    temp_in  = config%sea_surface_temperature
    salinity_in  = config%sea_surface_salinity

    write(*,'(a)') "Calculation Parameters:"
    write(*,'(a, a)') "  Mode: ", trim(config%mode)
    write(*,'(a, i0)') "  Number of angles:      ", n_angles
    write(*,'(a, i0)') "  Number of wind speeds: ", n_winds
    write(*,'(a, f6.1)') "  Sea surface temperature (K):       ", temp_in
    write(*,'(a, f6.1)') "  Sea surface salinity (PSU):       ", salinity_in

    if (trim(config%mode) == 'spectral') then
        ! For spectral mode, use first band settings
        ! Configuration automatically converts wavelength to wavenumber
        wavenum_start = config%wavenum_start(1)
        wavenum_end   = config%wavenum_end(1)
        n_wavenum     = config%n_wavenum(1)

        if (config%is_wavelength_input) then
            ! Wavelength input mode (already converted to wavenumber internally)
            write(*,'(a, f8.1)') "  Wavelength min: ", config%wavelen_min(1)
            write(*,'(a, f8.1)') "  Wavelength max: ", config%wavelen_max(1)
            write(*,'(a, i0)')   "  Number of wavelengths: ", n_wavenum
            write(*,'(a, f8.1, a, f8.1, a)') "  Converted to wavenumber range: ", &
                wavenum_start, " to ", wavenum_end, " cm⁻¹"
            write(*,'(a)') "  Using wavelength input mode"
        else
            ! Wavenumber input mode
            write(*,'(a, f8.1)') "  Wavenumber start: ", wavenum_start
            write(*,'(a, f8.1)') "  Wavenumber end:   ", wavenum_end
            write(*,'(a, i0)')   "  Number of wavenumbers: ", n_wavenum
            write(*,'(a)') "  Using wavenumber input mode"
        end if
    else if (trim(config%mode) == 'narrow-band') then
        write(*,'(a, i0)') "  Number of bands: ", config%n_bands
        do i = 1, config%n_bands
            write(*,'(a, i0, a, f8.1, a, f8.1, a, i0, a)') "    Band ", i, ": ", &
                config%wavenum_start(i), " to ", config%wavenum_end(i), &
                " cm-1 (", config%n_wavenum(i), " points)"
        end do
    else
        write(*,'(a, a)') "ERROR: Unknown mode: ", trim(config%mode)
        stop 1
    end if
    write(*,'(a)') ""

    !----------------------------------------------------------------------
    ! CROSS-PLATFORM OPTIMIZATION: Adaptive OpenMP configuration
    ! Optimized for both M-series (Apple Silicon) and Linux HPC systems
    !----------------------------------------------------------------------

    ! OpenMP environment variables should be set externally for best performance
    ! e.g., export OMP_NUM_THREADS = 14; export OMP_PROC_BIND = spread; export OMP_PLACES = cores
    write(*,'(a, i0, a)') "Using ", num_threads, &
        " OpenMP threads (set OMP_*variables externally for optimization)"

    !----------------------------------------------------------------------
    ! 3) Allocate arrays
    !----------------------------------------------------------------------
    if (trim(config%mode) == 'spectral') then
        ! For spectral mode, allocate arrays for the full spectral range
        allocate(wavenum(n_wavenum), wavelen(n_wavenum), real_n(n_wavenum), imag_n(n_wavenum), &
            refm(n_wavenum), angles(n_angles), winds(n_winds), &
            effective_angles(n_angles, n_winds), &
            emiss(n_angles, n_winds, n_wavenum), &
            emiss_noref(n_angles, n_winds, n_wavenum), &
            refl_cb(n_angles, n_winds, n_wavenum), stat = ierr)
    else
        ! For narrow-band mode, allocate for the largest band initially
        n_wavenum = maxval(config%n_wavenum(1:config%n_bands))
        allocate(wavenum(n_wavenum), real_n(n_wavenum), imag_n(n_wavenum), &
            refm(n_wavenum), angles(n_angles), winds(n_winds), &
            effective_angles(n_angles, n_winds), &
            emiss(n_angles, n_winds, n_wavenum), &
            emiss_noref(n_angles, n_winds, n_wavenum), &
            refl_cb(n_angles, n_winds, n_wavenum), stat = ierr)
    end if

    if (ierr /= 0) then
        print *, "Error allocating arrays"
        stop
    end if

    ! Initialize angles and winds from config
    do i = 1, n_angles
        angles(i) = config%view_angle(i)
    end do
    do i = 1, n_winds
        winds(i) = config%winds(i)
    end do

    !-------------------------------------------------------------------
    ! 4) Load effective and normal view angle LUTs, set up angle blending
    !-------------------------------------------------------------------
    if (config%use_effective_angle) then
        ! Angle blending strategy (Nalli et al. continuity correction):
        !   - Angles <= 70°: use effective_view_angle from effective_view_angle.nc
        !   - Angles >  70°: use normal view angles from normal_view_angle.nc
        !                     with continuity correction at 70° boundary
        !   - diff_70 = emiss_eff(70°) - emiss_nor(70°)
        !   - emiss_final(>70°) = emiss_nor(>70°) + diff_70
        !
        ! Read the effective_view_angle.nc (covers 0-70°)
        ncfile_name = 'data/effective_view_angle.nc'

        inquire(file = trim(ncfile_name), exist = file_exists)
        if (.not. file_exists) then
            print *, 'ERROR: NetCDF file not found: ', trim(ncfile_name)
            stop 1
        end if

        call file%open(trim(ncfile_name))
        if (.not. file%is_open()) then
            print *, 'ERROR: Failed to open NetCDF file: ', trim(ncfile_name)
            stop 1
        end if

        if (.not. file%exists('reference_angle')) then
            print *, 'ERROR: "reference_angle" not found in the effective_view_angle NETCDF file'
            call file%close()
            stop 1
        endif
        call file%get('reference_angle', reference_angle)

        if (.not. file%exists('reference_wind_speed')) then
            print *, 'ERROR: "reference_wind_speed" not found in the effective_view_angle NETCDF file'
            call file%close()
            stop 1
        endif
        call file%get('reference_wind_speed', reference_wind_speed)

        if (.not. file%exists('reference_sea_surface_temperature')) then
            print *, 'ERROR: "reference_sea_surface_temperature" not found in the effective_view_angle NETCDF file'
            call file%close()
            stop 1
        endif
        call file%get('reference_sea_surface_temperature', reference_sst)

        if (.not. file%exists('effective_view_angle')) then
            print *, 'ERROR: "effective_view_angle" not found in the effective_view_angle NETCDF file'
            call file%close()
            stop 1
        endif
        ! Read 3D LUT with ipermute=[3,2,1]: CDL (angle,wind,sst) → Fortran (sst,wind,angle) → permuted (angle,wind,sst)
        call file%get('effective_view_angle', effective_view_angle_lut_3d, ipermute=[3,2,1])
        call file%close()

        ! Interpolate 3D LUT to user SST → 2D slice (angle, wind)
        call interpolate_lut_sst(effective_view_angle_lut_3d, reference_sst, &
            config%sea_surface_temperature, effective_view_angle_lut)
        deallocate(effective_view_angle_lut_3d, reference_sst)

        ! Smooth the effective view angle LUT for C1 continuity
        ! Re-interpolates through 5-degree anchor points using monotone PCHIP
        call smooth_effective_angle_lut(reference_angle, reference_wind_speed, effective_view_angle_lut)

        ! Read the normal_view_angle.nc (covers 0-85°)
        ncfile_name = 'data/normal_view_angle.nc'

        inquire(file = trim(ncfile_name), exist = file_exists)
        if (.not. file_exists) then
            print *, 'ERROR: NetCDF file not found: ', trim(ncfile_name)
            stop 1
        end if

        call file%open(trim(ncfile_name))
        if (.not. file%is_open()) then
            print *, 'ERROR: Failed to open NetCDF file: ', trim(ncfile_name)
            stop 1
        end if

        if (.not. file%exists('reference_angle')) then
            print *, 'ERROR: "reference_angle" not found in the normal_view_angle NETCDF file'
            call file%close()
            stop 1
        endif
        call file%get('reference_angle', reference_angle_nor)

        if (.not. file%exists('reference_wind_speed')) then
            print *, 'ERROR: "reference_wind_speed" not found in the normal_view_angle NETCDF file'
            call file%close()
            stop 1
        endif
        call file%get('reference_wind_speed', reference_wind_speed_nor)

        if (.not. file%exists('effective_view_angle')) then
            print *, 'ERROR: "effective_view_angle" not found in the normal_view_angle NETCDF file'
            call file%close()
            stop 1
        endif
        call file%get('effective_view_angle', normal_view_angle_lut)
        call file%close()

        !-------------------------------------------------------------------
        ! Determine angle splitting at 70° boundary
        !-------------------------------------------------------------------
        n_angles_le70 = 0
        n_angles_gt70 = 0
        has_70_exact = .false.
        do i = 1, n_angles
            if (angles(i) <= 70.0_jprd + 1.0e-6_jprd) then
                n_angles_le70 = n_angles_le70 + 1
                if (abs(angles(i) - 70.0_jprd) < 1.0e-6_jprd) has_70_exact = .true.
            else
                n_angles_gt70 = n_angles_gt70 + 1
            end if
        end do
        has_angles_gt70 = (n_angles_gt70 > 0)

        if (has_angles_gt70) then
            ! Build angles_eff_calc: user angles <= 70° plus 70° bridge if not present
            if (has_70_exact) then
                n_eff_calc = n_angles_le70
            else
                n_eff_calc = n_angles_le70 + 1  ! Add 70° as bridge point
            end if
            allocate(angles_eff_calc(n_eff_calc))
            idx_70_eff = 0
            do i = 1, n_angles
                if (angles(i) <= 70.0_jprd + 1.0e-6_jprd) then
                    idx_70_eff = idx_70_eff + 1
                    angles_eff_calc(idx_70_eff) = angles(i)
                end if
            end do
            if (.not. has_70_exact) then
                angles_eff_calc(n_eff_calc) = 70.0_jprd
            end if
            idx_70_eff = n_eff_calc  ! 70° is always last in this array

            ! Build angles_nor_calc: 70° bridge + user angles > 70°
            n_nor_calc = 1 + n_angles_gt70
            allocate(angles_nor_calc(n_nor_calc))
            angles_nor_calc(1) = 70.0_jprd
            k = 1
            do i = 1, n_angles
                if (angles(i) > 70.0_jprd + 1.0e-6_jprd) then
                    k = k + 1
                    angles_nor_calc(k) = angles(i)
                end if
            end do

            ! Compute effective angles for each group
            allocate(eff_angles_eff(n_eff_calc, n_winds))
            allocate(eff_angles_nor(n_nor_calc, n_winds))

            call effective_view_angle(reference_angle, reference_wind_speed, &
                effective_view_angle_lut, angles_eff_calc, winds, eff_angles_eff)

            call effective_view_angle(reference_angle_nor, reference_wind_speed_nor, &
                normal_view_angle_lut, angles_nor_calc, winds, eff_angles_nor)
        else
            ! All angles <= 70°: use effective_view_angle only (original behavior)
            call effective_view_angle(reference_angle, reference_wind_speed, &
                effective_view_angle_lut, angles, winds, effective_angles)
        end if

    else
        ! use_effective_angle = .false.: use only normal view angle LUT for all angles
        write(*,'(a)') "Using normal view angles (effective angle LUT disabled)"

        ncfile_name = 'data/normal_view_angle.nc'

        inquire(file = trim(ncfile_name), exist = file_exists)
        if (.not. file_exists) then
            print *, 'ERROR: NetCDF file not found: ', trim(ncfile_name)
            stop 1
        end if

        call file%open(trim(ncfile_name))
        if (.not. file%is_open()) then
            print *, 'ERROR: Failed to open NetCDF file: ', trim(ncfile_name)
            stop 1
        end if

        if (.not. file%exists('reference_angle')) then
            print *, 'ERROR: "reference_angle" not found in the normal_view_angle NETCDF file'
            call file%close()
            stop 1
        endif
        call file%get('reference_angle', reference_angle_nor)

        if (.not. file%exists('reference_wind_speed')) then
            print *, 'ERROR: "reference_wind_speed" not found in the normal_view_angle NETCDF file'
            call file%close()
            stop 1
        endif
        call file%get('reference_wind_speed', reference_wind_speed_nor)

        if (.not. file%exists('effective_view_angle')) then
            print *, 'ERROR: "effective_view_angle" not found in the normal_view_angle NETCDF file'
            call file%close()
            stop 1
        endif
        call file%get('effective_view_angle', normal_view_angle_lut)
        call file%close()

        ! No angle blending — all angles use normal LUT
        has_angles_gt70 = .false.

        call effective_view_angle(reference_angle_nor, reference_wind_speed_nor, &
            normal_view_angle_lut, angles, winds, effective_angles)
    end if

    print *, ""
    print *, "Input configurations:"
    write(*,'(A)', advance='no') "  angles (degrees): "
    do i = 1, n_angles
        write(*,'(F6.1)', advance='no') angles(i)
        if (i < n_angles) write(*,'(A)', advance='no') ', '
    end do
    write(*,*)
    write(*,'(A)', advance='no') "  winds (m/s): "
    do i = 1, n_winds
        write(*,'(F5.1)', advance='no') winds(i)
        if (i < n_winds) write(*,'(A)', advance='no') ', '
    end do
    write(*,*)
    print *, ""

    !----------------------------------------------------------------------
    ! 5) Now branch into spectral vs narrow-band processing with polarization support
    !----------------------------------------------------------------------
    write(*,'(a, a)') "  Polarization mode: ", trim(config%polarization)
    
    if (trim(config%mode) == 'spectral') then
        if (trim(config%polarization) == 'enable') then
            ! POLARIZED SPECTRAL MODE: Calculate vertical and horizontal polarized emissivities
            call process_polarized_spectral_mode()
        else
            ! SPECTRAL MODE: Calculate unpolarized emissivity for full spectral range
            call process_spectral_mode()
        end if
    else if (trim(config%mode) == 'narrow-band') then
        if (trim(config%polarization) == 'enable') then
            write(*,'(a)') "ERROR: Polarized narrow-band mode not yet implemented"
            stop 1
        else
            ! NARROW-BAND MODE: Calculate for each band and integrate with Planck function
            call process_narrow_band_mode()
        end if
    else
        write(*,'(a, a)') "ERROR: Unsupported mode: ", trim(config%mode)
        stop 1
    end if

contains

    !----------------------------------------------------------------------
    ! Helper subroutine: to_lowercase
    ! Converts a string to lowercase for case-insensitive comparison
    !----------------------------------------------------------------------
    subroutine to_lowercase(str)
        character(len=*), intent(inout) :: str
        integer :: i, ic
        do i = 1, len_trim(str)
            ic = ichar(str(i:i))
            if (ic >= 65 .and. ic <= 90) str(i:i) = char(ic + 32)
        end do
    end subroutine to_lowercase

    !----------------------------------------------------------------------
    ! Process spectral mode (with angle blending at 70° boundary)
    !----------------------------------------------------------------------
    subroutine process_spectral_mode()
        ! Local variables for computation
        integer(kind = jpim):: idx
        logical:: ncfile_out_exists
        real(kind = jprd), allocatable:: emiss_noref_local(:,:,:)  ! For no_multiple_reflection mode
        ! Local variables for angle blending
        real(kind = jprd), allocatable:: emiss_eff_part(:,:,:)     ! Emissivity from effective angles
        real(kind = jprd), allocatable:: emiss_nor_part(:,:,:)     ! Emissivity from normal angles
        real(kind = jprd), allocatable:: emiss_noref_eff(:,:,:)    ! No-refl emissivity from eff angles
        real(kind = jprd), allocatable:: emiss_noref_nor(:,:,:)    ! No-refl emissivity from normal angles

        ! Get refractive indices (complex) from NetCDF file
        call get_sea_water_nk(config, wavenum_start, wavenum_end, n_wavenum, temp_in, salinity_in, &
            real_n, imag_n, wavenum)

        ! If wavelength input, generate wavelength array from wavenumber for output
        if (config%is_wavelength_input) then
            do i = 1, n_wavenum
                wavelen(i) = 10000.0_jprd/wavenum(i)
            end do
        end if

        refm = cmplx(real_n, imag_n, kind = jprd)
        write(*,'(a, a)') "Complex refractive index loaded from: ", trim(config%refractive_index_filename)
        if (n_wavenum > 0) then
            if (config%is_wavelength_input) then
                write(*,'(a, f8.2, a, f8.2, a)') "Wavelength range: ", &
                    wavelen(1), " to ", wavelen(n_wavenum), " μm"
                write(*,'(a, f8.2, a, f8.2, a)') "Converted wavenumber range: ", &
                    wavenum(1), " to ", wavenum(n_wavenum), " cm⁻¹"
            else
                write(*,'(a, f8.2, a, f8.2, a)') "Wavenumber range: ", &
                    wavenum(1), " to ", wavenum(n_wavenum), " cm⁻¹"
            end if
        else
            write(*,'(a)') "Spectral range: [not available-no data read]"
        end if
        write(*,'(a)') ""

        ! Perform emissivity calculations
        call system_clock(start_clock, clock_rate)

        write(*,'(a)') "Starting exact accuracy emissivity calculations..."
        write(*,'(a, i0, a, i0, a, i0, a, i0)') "Total computations: ", n_wavenum, " x ", n_angles, &
            " x ", n_winds, " = ", n_wavenum*n_angles*n_winds
        write(*,'(a)') ""

        if (has_angles_gt70) then
            !----------------------------------------------------------
            ! Angle blending: dual computation + continuity correction
            !----------------------------------------------------------
            allocate(emiss_eff_part(n_eff_calc, n_winds, n_wavenum))
            allocate(emiss_nor_part(n_nor_calc, n_winds, n_wavenum))

            if (config%no_multiple_reflection) then
                allocate(emiss_noref_local(n_angles, n_winds, n_wavenum))
                allocate(emiss_noref_eff(n_eff_calc, n_winds, n_wavenum))
                allocate(emiss_noref_nor(n_nor_calc, n_winds, n_wavenum))
                call get_emissivity_optimized(refm, eff_angles_eff, winds, emiss_eff_part, emiss_noref_eff)
                call get_emissivity_optimized(refm, eff_angles_nor, winds, emiss_nor_part, emiss_noref_nor)
                call blend_angle_results(emiss_noref_eff, emiss_noref_nor, emiss_noref_local, n_wavenum)
                deallocate(emiss_noref_eff, emiss_noref_nor)
            else
                call get_emissivity_optimized(refm, eff_angles_eff, winds, emiss_eff_part)
                call get_emissivity_optimized(refm, eff_angles_nor, winds, emiss_nor_part)
            end if

            call blend_angle_results(emiss_eff_part, emiss_nor_part, emiss, n_wavenum)
            deallocate(emiss_eff_part, emiss_nor_part)
        else
            ! All angles <= 70°: single computation (original behavior)
            if (config%no_multiple_reflection) then
                allocate(emiss_noref_local(n_angles, n_winds, n_wavenum))
                call get_emissivity_optimized(refm, effective_angles, winds, emiss, emiss_noref_local)
            else
                call get_emissivity_optimized(refm, effective_angles, winds, emiss)
            end if
        end if

        ! Smooth slope discontinuity at 70° boundary, then enforce monotonicity
        if (has_angles_gt70) call smooth_70_boundary(emiss)
        call enforce_emissivity_monotonicity(emiss)
        if (config%no_multiple_reflection) then
            if (has_angles_gt70) call smooth_70_boundary(emiss_noref_local)
            call enforce_emissivity_monotonicity(emiss_noref_local)
        end if

        call system_clock(end_clock)
        elapsed_time = real(end_clock-start_clock, jprd) / real(clock_rate, jprd)
        performance_rate = real(n_wavenum*n_angles*n_winds, jprd) / elapsed_time

        write(*,'(a)') "Calculation completed successfully!"
        write(*,'(a, f8.2, a)') "  Total time: ", elapsed_time, " seconds"
        write(*,'(a, f12.0, a)') "  Performance: ", performance_rate, " calculations/second"
        write(*,'(a)') ""

        ! Check if output file exists and remove it if necessary
        inquire(file = trim(output_file), exist = ncfile_out_exists)
        if (ncfile_out_exists) then
            write(*,'(a, a)') "|----- Deleting existing output file: ", trim(output_file)
            call execute_command_line("rm -f " // trim(output_file))
        end if

        ! Save spectral results
        ! Save appropriate emissivity based on no_multiple_reflection config
        write(*,'(a, a)') "Writing results to NetCDF file: ", trim(output_file)
        if (config%no_multiple_reflection) then
            write(*,'(a)') "  Output mode: emissivity WITHOUT multiple reflection (first-order Fresnel)"
            if (config%is_wavelength_input) then
                call save_spectral_emissivity_wl(output_file, config, angles, winds, wavelen, &
                    emiss_noref_local, temp_in, salinity_in)
            else
                call save_spectral_emissivity(output_file, config, angles, winds, wavenum, &
                    emiss_noref_local, temp_in, salinity_in)
            end if
            deallocate(emiss_noref_local)
        else
            if (config%is_wavelength_input) then
                call save_spectral_emissivity_wl(output_file, config, angles, winds, wavelen, &
                    emiss, temp_in, salinity_in)
            else
                call save_spectral_emissivity(output_file, config, angles, winds, wavenum, &
                    emiss, temp_in, salinity_in)
            end if
        end if

    end subroutine process_spectral_mode

    !----------------------------------------------------------------------
    ! Process polarized spectral mode (V and H polarizations)
    !----------------------------------------------------------------------
    subroutine process_polarized_spectral_mode()
        ! Local variables for polarized computation
        real(kind = jprd), allocatable:: emiss_v(:,:,:), emiss_h(:,:,:)
        real(kind = jprd), allocatable:: emiss_noref_v(:,:,:), emiss_noref_h(:,:,:)
        real(kind = jprd), allocatable:: refl_cb_v(:,:,:), refl_cb_h(:,:,:)
        ! Local variables for angle blending (polarized)
        real(kind = jprd), allocatable:: emiss_v_eff(:,:,:), emiss_h_eff(:,:,:)
        real(kind = jprd), allocatable:: emiss_v_nor(:,:,:), emiss_h_nor(:,:,:)
        real(kind = jprd), allocatable:: noref_v_eff(:,:,:), noref_h_eff(:,:,:)
        real(kind = jprd), allocatable:: noref_v_nor(:,:,:), noref_h_nor(:,:,:)
        real(kind = jprd), allocatable:: refl_v_eff(:,:,:), refl_h_eff(:,:,:)
        real(kind = jprd), allocatable:: refl_v_nor(:,:,:), refl_h_nor(:,:,:)
        integer(kind = jpim):: idx
        logical:: ncfile_out_exists
        real(kind=jprd) :: elapsed
        integer(kind=jpim) :: start_clock, end_clock, clock_rate

        ! Get refractive indices (complex) from NetCDF file
        call get_sea_water_nk(config, wavenum_start, wavenum_end, n_wavenum, temp_in, salinity_in, &
            real_n, imag_n, wavenum)

        ! If wavelength input, generate wavelength array from wavenumber for output
        if (config%is_wavelength_input) then
            do i = 1, n_wavenum
                wavelen(i) = 10000.0_jprd/wavenum(i)
            end do
        end if

        refm = cmplx(real_n, imag_n, kind = jprd)
        write(*,'(a, a)') "Complex refractive index loaded from: ", trim(config%refractive_index_filename)
        if (n_wavenum > 0) then
            if (config%is_wavelength_input) then
                write(*,'(a, f8.2, a, f8.2, a)') "Wavelength range: ", &
                    wavelen(1), " to ", wavelen(n_wavenum), " μm"
                write(*,'(a, f8.2, a, f8.2, a)') "Converted wavenumber range: ", &
                    wavenum(1), " to ", wavenum(n_wavenum), " cm⁻¹"
            else
                write(*,'(a, f8.2, a, f8.2, a)') "Wavenumber range: ", &
                    wavenum(1), " to ", wavenum(n_wavenum), " cm⁻¹"
            end if
        else
            write(*,'(a)') "Spectral range: [not available-no data read]"
        end if
        write(*,'(a)') ""

        ! Perform polarized emissivity calculations
        call system_clock(start_clock, clock_rate)
        write(*,'(a)') "Starting polarized emissivity calculations (V and H components)..."
        write(*,'(a, i0, a, i0, a, i0, a, i0)') "Total computations: ", n_wavenum, " x ", n_angles, &
            " x ", n_winds, " = ", n_wavenum*n_angles*n_winds
        write(*,'(a)') ""

        ! Pre-allocate final output arrays
        allocate(emiss_v(n_angles, n_winds, n_wavenum), &
                 emiss_h(n_angles, n_winds, n_wavenum), &
                 emiss_noref_v(n_angles, n_winds, n_wavenum), &
                 emiss_noref_h(n_angles, n_winds, n_wavenum), &
                 refl_cb_v(n_angles, n_winds, n_wavenum), &
                 refl_cb_h(n_angles, n_winds, n_wavenum))

        if (has_angles_gt70) then
            !----------------------------------------------------------
            ! Angle blending for polarized mode
            !----------------------------------------------------------
            allocate(emiss_v_eff(n_eff_calc, n_winds, n_wavenum), &
                     emiss_h_eff(n_eff_calc, n_winds, n_wavenum), &
                     noref_v_eff(n_eff_calc, n_winds, n_wavenum), &
                     noref_h_eff(n_eff_calc, n_winds, n_wavenum), &
                     refl_v_eff(n_eff_calc, n_winds, n_wavenum), &
                     refl_h_eff(n_eff_calc, n_winds, n_wavenum))
            allocate(emiss_v_nor(n_nor_calc, n_winds, n_wavenum), &
                     emiss_h_nor(n_nor_calc, n_winds, n_wavenum), &
                     noref_v_nor(n_nor_calc, n_winds, n_wavenum), &
                     noref_h_nor(n_nor_calc, n_winds, n_wavenum), &
                     refl_v_nor(n_nor_calc, n_winds, n_wavenum), &
                     refl_h_nor(n_nor_calc, n_winds, n_wavenum))

            call get_polarized_emissivity_optimized(refm, eff_angles_eff, winds, &
                emiss_v_eff, noref_v_eff, refl_v_eff, emiss_h_eff, noref_h_eff, refl_h_eff)
            call get_polarized_emissivity_optimized(refm, eff_angles_nor, winds, &
                emiss_v_nor, noref_v_nor, refl_v_nor, emiss_h_nor, noref_h_nor, refl_h_nor)

            ! Blend each polarization component
            call blend_angle_results(emiss_v_eff, emiss_v_nor, emiss_v, n_wavenum)
            call blend_angle_results(emiss_h_eff, emiss_h_nor, emiss_h, n_wavenum)
            call blend_angle_results(noref_v_eff, noref_v_nor, emiss_noref_v, n_wavenum)
            call blend_angle_results(noref_h_eff, noref_h_nor, emiss_noref_h, n_wavenum)
            call blend_angle_results(refl_v_eff, refl_v_nor, refl_cb_v, n_wavenum)
            call blend_angle_results(refl_h_eff, refl_h_nor, refl_cb_h, n_wavenum)

            deallocate(emiss_v_eff, emiss_h_eff, noref_v_eff, noref_h_eff, refl_v_eff, refl_h_eff)
            deallocate(emiss_v_nor, emiss_h_nor, noref_v_nor, noref_h_nor, refl_v_nor, refl_h_nor)
        else
            ! All angles <= 70°: single computation (original behavior)
            call get_polarized_emissivity_optimized(refm, effective_angles, winds, &
                emiss_v, emiss_noref_v, refl_cb_v, emiss_h, emiss_noref_h, refl_cb_h)
        end if

        ! Smooth slope discontinuity at 70° boundary, then enforce monotonicity
        if (has_angles_gt70) then
            call smooth_70_boundary(emiss_v)
            call smooth_70_boundary(emiss_h)
        end if
        ! Note: Do NOT enforce monotonicity on V-pol (emiss_v).
        ! V-pol emissivity physically increases from nadir toward Brewster angle
        ! (~53° for water) before decreasing at large angles. Clamping destroys this.
        ! H-pol emissivity is physically monotonically decreasing with angle.
        call enforce_emissivity_monotonicity(emiss_h)

        call system_clock(end_clock)
        elapsed = real(end_clock - start_clock, jprd) / real(clock_rate, jprd)

        write(*,'(a)') ""
        write(*,'(a)') "=== POLARIZED EMISSIVITY CALCULATION COMPLETED ==="
        write(*,'(a, f8.3, a)') "Elapsed time: ", elapsed, " seconds"
        write(*,'(a, f8.1, a)') "Performance:  ", real(n_angles*n_winds*n_wavenum, jprd)/elapsed, " calculations/second"
        write(*,'(a)') ""

        ! Check for output file existence and save results
        inquire(file=trim(output_file), exist=ncfile_out_exists)
        if (ncfile_out_exists) then
            write(*,'(a, a)') "Warning: Output file already exists and will be overwritten: ", trim(output_file)
        end if

        ! Save polarized results
        write(*,'(a, a)') "Writing polarized results to NetCDF file: ", trim(output_file)
        if (config%is_wavelength_input) then
            call save_spectral_polarized_emissivity_wl(output_file, config, angles, winds, wavelen, &
                emiss_v, emiss_h, temp_in, salinity_in)
        else
            call save_spectral_polarized_emissivity(output_file, config, angles, winds, wavenum, &
                emiss_v, emiss_h, temp_in, salinity_in)
        end if
        write(*,'(a)') "Polarized emissivity computation and saving completed successfully."

        deallocate(emiss_v, emiss_h, emiss_noref_v, emiss_noref_h, refl_cb_v, refl_cb_h)
    end subroutine process_polarized_spectral_mode


    !----------------------------------------------------------------------
    ! Process narrow-band mode with Planck integration
    !
    ! Major Performance Optimizations:
    ! 1. Single NetCDF read for entire spectrum (eliminates 15 redundant I/O ops)
    ! 2. Single get_emissivity_optimized call for ALL winds (eliminates serial loop)
    ! 3. Uses unified wind processing architecture matching polarized mode
    !----------------------------------------------------------------------
    subroutine process_narrow_band_mode()
        use mathlib, only: planck_radiance, get_trapz, error_type, integrate_planck_vectorized

        ! Local variables for narrow-band processing
        real(kind = jprd), allocatable:: emiss_narrow(:,:,:)
        real(kind = jprd), allocatable:: wavenum_band(:)

        ! Full spectrum arrays for single NetCDF read
        real(kind = jprd), allocatable:: full_real_n(:), full_imag_n(:), full_wavenum(:)
        real(kind = jprd), allocatable:: full_emiss(:,:,:)  ! (n_angles, n_winds, full_n_points)
        real(kind = jprd), allocatable:: full_emiss_noref(:,:,:)
        complex(kind = jprd), allocatable:: full_refm(:)
        real(kind = jprd):: full_wavenum_start, full_wavenum_end
        integer(kind = jpim):: full_n_points, start_idx

        ! Local variables for angle blending in narrow-band mode
        real(kind = jprd), allocatable:: full_emiss_eff(:,:,:)
        real(kind = jprd), allocatable:: full_emiss_nor(:,:,:)
        real(kind = jprd), allocatable:: full_emiss_noref_eff(:,:,:)
        real(kind = jprd), allocatable:: full_emiss_noref_nor(:,:,:)

        real(kind = jprd):: band_start
        integer(kind = jpim):: band_idx, n_band_points
        type(error_type):: err_trapz
        logical:: ncfile_out_exists

        ! Allocate arrays for narrow-band results
        allocate(emiss_narrow(n_angles, n_winds, config%n_bands))

        write(*,'(a)') "Starting narrow-band emissivity calculations..."
        write(*,'(a, i0, a)') "Processing ", config%n_bands, " narrow bands"
        write(*,'(a)') ""

        ! Calculate full spectrum range across all bands
        full_wavenum_start = minval(config%wavenum_start(1:config%n_bands))
        full_wavenum_end = maxval(config%wavenum_end(1:config%n_bands))
        full_n_points = sum(config%n_wavenum(1:config%n_bands)) - (config%n_bands-1)

        write(*,'(a)') ""
        write(*,'(a, f8.1, a, f8.1, a, i0, a)') "Full spectrum: ", full_wavenum_start, " to ", &
            full_wavenum_end, " cm-1 (", full_n_points, " total points)"
        write(*,'(a)') ""

        ! Read full spectrum data once
        call get_sea_water_nk(config, full_wavenum_start, full_wavenum_end, full_n_points, &
            temp_in, salinity_in, full_real_n, full_imag_n, full_wavenum)

        allocate(full_refm(full_n_points))
        full_refm = cmplx(full_real_n, full_imag_n, kind = jprd)

        allocate(wavenum_band(maxval(config%n_wavenum(1:config%n_bands))))

        call system_clock(start_clock, clock_rate)

        write(*,'(a, i0, a, i0, a, i0, a, i0)') "  Total computations: ", n_angles, " x ", n_winds, &
            " x ", full_n_points, " = ", n_angles * n_winds * full_n_points

        ! Pre-allocate full emissivity result array
        allocate(full_emiss(n_angles, n_winds, full_n_points))

        if (has_angles_gt70) then
            !----------------------------------------------------------
            ! Angle blending for narrow-band mode
            !----------------------------------------------------------
            allocate(full_emiss_eff(n_eff_calc, n_winds, full_n_points))
            allocate(full_emiss_nor(n_nor_calc, n_winds, full_n_points))

            if (config%no_multiple_reflection) then
                allocate(full_emiss_noref(n_angles, n_winds, full_n_points))
                allocate(full_emiss_noref_eff(n_eff_calc, n_winds, full_n_points))
                allocate(full_emiss_noref_nor(n_nor_calc, n_winds, full_n_points))
                call get_emissivity_optimized(full_refm, eff_angles_eff, winds, full_emiss_eff, full_emiss_noref_eff)
                call get_emissivity_optimized(full_refm, eff_angles_nor, winds, full_emiss_nor, full_emiss_noref_nor)
                call blend_angle_results(full_emiss_noref_eff, full_emiss_noref_nor, full_emiss_noref, full_n_points)
                deallocate(full_emiss_noref_eff, full_emiss_noref_nor)
            else
                call get_emissivity_optimized(full_refm, eff_angles_eff, winds, full_emiss_eff)
                call get_emissivity_optimized(full_refm, eff_angles_nor, winds, full_emiss_nor)
            end if

            call blend_angle_results(full_emiss_eff, full_emiss_nor, full_emiss, full_n_points)
            deallocate(full_emiss_eff, full_emiss_nor)
        else
            ! All angles <= 70°: single computation (original behavior)
            if (config%no_multiple_reflection) then
                allocate(full_emiss_noref(n_angles, n_winds, full_n_points))
                call get_emissivity_optimized(full_refm, effective_angles, winds, full_emiss, full_emiss_noref)
            else
                call get_emissivity_optimized(full_refm, effective_angles, winds, full_emiss)
            end if
        end if

        ! Smooth slope discontinuity at 70° boundary, then enforce monotonicity
        if (has_angles_gt70) call smooth_70_boundary(full_emiss)
        call enforce_emissivity_monotonicity(full_emiss)
        if (config%no_multiple_reflection .and. allocated(full_emiss_noref)) then
            if (has_angles_gt70) call smooth_70_boundary(full_emiss_noref)
            call enforce_emissivity_monotonicity(full_emiss_noref)
        end if

        ! Extract results for each band from full spectrum and perform Planck integration
        do band_idx = 1, config%n_bands
            band_start = config%wavenum_start(band_idx)
            n_band_points = config%n_wavenum(band_idx)

            ! Find start index in full spectrum for this band
            start_idx = 1
            do i = 1, full_n_points
                if (abs(full_wavenum(i) - band_start) < abs(full_wavenum(start_idx) - band_start)) then
                    start_idx = i
                end if
            end do

            wavenum_band(1:n_band_points) = full_wavenum(start_idx:start_idx+n_band_points-1)

            do j = 1, n_winds
                if (config%no_multiple_reflection) then
                    call integrate_planck_vectorized(wavenum_band(1:n_band_points), &
                        full_emiss_noref(:, j:j, start_idx:start_idx+n_band_points-1), temp_in, &
                        n_angles, 1, n_band_points, emiss_narrow(:, j:j, band_idx), err_trapz)
                else
                    call integrate_planck_vectorized(wavenum_band(1:n_band_points), &
                        full_emiss(:, j:j, start_idx:start_idx+n_band_points-1), temp_in, &
                        n_angles, 1, n_band_points, emiss_narrow(:, j:j, band_idx), err_trapz)
                end if
            end do
        end do

        call system_clock(end_clock)
        elapsed_time = real(end_clock-start_clock, jprd) / real(clock_rate, jprd)

        write(*,'(a)') "Narrow-band calculations completed!"
        write(*,'(a, f8.2, a)') "  Total time: ", elapsed_time, " seconds"
        write(*,'(a)') ""

        ! Check if output file exists and remove it if necessary
        inquire(file = trim(output_file), exist = ncfile_out_exists)
        if (ncfile_out_exists) then
            write(*,'(a, a)') "|----- Deleting existing output file: ", trim(output_file)
            call execute_command_line("rm -f " // trim(output_file))
        end if

        write(*,'(a, a)') "Writing narrow-band results to NetCDF file: ", trim(output_file)
        if (config%no_multiple_reflection) then
            write(*,'(a)') "  Output mode: emissivity WITHOUT multiple reflection (first-order Fresnel)"
        end if
        call save_narrow_band_emissivity(output_file, config, angles, winds, &
            config%wavenum_start(1:config%n_bands), &
            config%wavenum_end(1:config%n_bands), emiss_narrow, &
            temp_in, salinity_in)

        ! Clean up arrays
        deallocate(emiss_narrow)
        if (allocated(full_real_n)) deallocate(full_real_n)
        if (allocated(full_imag_n)) deallocate(full_imag_n)
        if (allocated(full_wavenum)) deallocate(full_wavenum)
        if (allocated(full_refm)) deallocate(full_refm)
        if (allocated(full_emiss)) deallocate(full_emiss)
        if (allocated(full_emiss_noref)) deallocate(full_emiss_noref)
        if (allocated(wavenum_band)) deallocate(wavenum_band)

    end subroutine process_narrow_band_mode

    !----------------------------------------------------------------------
    ! blend_angle_results: Merge effective-angle and normal-angle emissivities
    ! with continuity correction at 70° boundary
    !
    ! emiss_eff_part: emissivity computed with effective angles (n_eff_calc, n_winds, n_spec)
    !                 70° is at index idx_70_eff (last element)
    ! emiss_nor_part: emissivity computed with normal angles (n_nor_calc, n_winds, n_spec)
    !                 70° is at index 1 (first element)
    ! emiss_out:      blended result (n_angles, n_winds, n_spec)
    ! n_spec:         number of spectral points (wavenumbers or bands)
    !----------------------------------------------------------------------
    subroutine blend_angle_results(emiss_eff_part, emiss_nor_part, emiss_out, n_spec)
        real(kind = jprd), intent(in)  :: emiss_eff_part(:,:,:)
        real(kind = jprd), intent(in)  :: emiss_nor_part(:,:,:)
        real(kind = jprd), intent(out) :: emiss_out(:,:,:)
        integer(kind = jpim), intent(in) :: n_spec

        real(kind = jprd), allocatable :: diff_70(:,:)  ! (n_winds, n_spec)
        integer(kind = jpim) :: ii, jj, kk, ie, in

        ! Compute continuity correction: diff_70 = emiss_eff(70°) - emiss_nor(70°)
        allocate(diff_70(n_winds, n_spec))
        do kk = 1, n_spec
            do jj = 1, n_winds
                diff_70(jj, kk) = emiss_eff_part(idx_70_eff, jj, kk) &
                    - emiss_nor_part(1, jj, kk)
            end do
        end do

        ! Map results to output array
        ie = 0   ! Index into emiss_eff_part (angles <= 70°)
        in = 1   ! Index into emiss_nor_part (skip 70° bridge at index 1)
        do ii = 1, n_angles
            if (angles(ii) <= 70.0_jprd + 1.0e-6_jprd) then
                ! Angles <= 70°: use effective-angle emissivity directly
                ie = ie + 1
                emiss_out(ii, :, :) = emiss_eff_part(ie, :, :)
            else
                ! Angles > 70°: use normal-angle emissivity + continuity correction
                in = in + 1
                do kk = 1, n_spec
                    do jj = 1, n_winds
                        emiss_out(ii, jj, kk) = emiss_nor_part(in, jj, kk) &
                            + diff_70(jj, kk)
                    end do
                end do
            end if
        end do

        deallocate(diff_70)
    end subroutine blend_angle_results


    !----------------------------------------------------------------------
    ! Smooth emissivity around the 70° model boundary to eliminate slope
    ! discontinuity between effective-angle and normal-angle regimes.
    ! Uses multiple passes of 3-point weighted averaging in [66°, 74°].
    !----------------------------------------------------------------------
    subroutine smooth_70_boundary(emiss_array)
        real(kind=jprd), intent(inout) :: emiss_array(:,:,:)  ! (n_angles, n_winds, n_spec)
        integer(kind=jpim) :: na, nw, ns, i, j, k, pass
        real(kind=jprd), parameter :: SMOOTH_LO = 66.0_jprd
        real(kind=jprd), parameter :: SMOOTH_HI = 74.0_jprd
        integer(kind=jpim), parameter :: N_PASSES = 5
        real(kind=jprd), allocatable :: temp_col(:)

        na = size(emiss_array, 1)
        nw = size(emiss_array, 2)
        ns = size(emiss_array, 3)

        allocate(temp_col(na))

        do k = 1, ns
            do j = 1, nw
                do pass = 1, N_PASSES
                    temp_col(:) = emiss_array(:,j,k)
                    do i = 2, na - 1
                        if (angles(i) >= SMOOTH_LO .and. angles(i) <= SMOOTH_HI) then
                            emiss_array(i,j,k) = 0.25_jprd * temp_col(i-1) + &
                                0.50_jprd * temp_col(i) + &
                                0.25_jprd * temp_col(i+1)
                        end if
                    end do
                end do
            end do
        end do

        deallocate(temp_col)
    end subroutine smooth_70_boundary


    !----------------------------------------------------------------------
    ! Enforce monotonic decrease of emissivity with increasing view angle.
    ! Forward-clamps: if emiss(i) > emiss(i-1), set emiss(i) = emiss(i-1).
    ! This is a lightweight safety net for any residual non-monotonicity.
    !----------------------------------------------------------------------
    subroutine enforce_emissivity_monotonicity(emiss_array)
        real(kind=jprd), intent(inout) :: emiss_array(:,:,:)  ! (n_angles, n_winds, n_spec)
        integer(kind=jpim) :: na, nw, ns, i, j, k, n_clamped

        na = size(emiss_array, 1)
        nw = size(emiss_array, 2)
        ns = size(emiss_array, 3)

        n_clamped = 0
        do k = 1, ns
            do j = 1, nw
                do i = 2, na
                    if (emiss_array(i,j,k) > emiss_array(i-1,j,k)) then
                        emiss_array(i,j,k) = emiss_array(i-1,j,k)
                        n_clamped = n_clamped + 1
                    end if
                end do
            end do
        end do
        if (n_clamped > 0) then
            write(*,'(a,i0,a)') "  Monotonicity safety net: clamped ", n_clamped, " values"
        end if
    end subroutine enforce_emissivity_monotonicity

end program driver
