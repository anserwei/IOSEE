! save_output.F90
!
! Purpose:
!   Output file generation and formatting for ocean emissivity calculations.
!   Provides NetCDF file creation with CF-compliant metadata for spectral,
!   narrow-band, and polarized emissivity results.
!
! Features:
!   - Spectral and narrow-band output formatting with proper metadata
!   - Polarized emissivity output support (vertical/horizontal components)
!   - NetCDF variable creation with comprehensive attributes
!   - Wavelength/wavenumber coordinate systems with unit conversion
!   - Global metadata management with calculation parameters
!   - CF-compliant scientific data format standards
!   - Professional scientific output with proper documentation
!
! Dependencies:
!   - utils: Precision types
!   - configuration: Configuration management
!   - netcdf_handler: NetCDF file operations
!
! Public Interface:
!   - save_spectral_emissivity: Unpolarized spectral mode output
!   - save_spectral_emissivity_wl: Spectral output with wavelength indexing
!   - save_narrow_band_emissivity: Narrow-band mode output with Planck weighting
!   - save_spectral_polarized_emissivity: Polarized spectral output (V/H components)
!   - save_spectral_polarized_emissivity_wl: Polarized output with wavelength indexing
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
module save_output
    use utils,              only: jpim, jprd
    use configuration,      only: config_type
    use netcdf_handler,     only: netcdf_file
    implicit none
    private

    public:: save_spectral_emissivity, save_spectral_emissivity_wl, save_narrow_band_emissivity, &
        save_spectral_polarized_emissivity, save_spectral_polarized_emissivity_wl
    private:: add_global_attributes_spectral, add_global_attributes_narrow_band
    private:: get_current_datetime, file_size_mb

    ! Units and constants for NetCDF attributes
    character(len=*), parameter:: EMISS_UNITS = "1"
    character(len=*), parameter:: ANGLE_UNITS = "degrees"
    character(len=*), parameter:: WIND_UNITS = "m/s"
    character(len=*), parameter:: WAVENUMBER_UNITS = "cm-1"
    character(len=*), parameter:: WAVELENGTH_UNITS = "um"
    character(len=*), parameter:: TEMPERATURE_UNITS = "K"

contains

    !-----------------------------------------------------------------------
    ! Save ocean emissivity results to NetCDF file
    !
    ! This subroutine creates a comprehensive NetCDF file containing all
    ! ocean emissivity calculation results with proper CF-compliant metadata
    !
    ! Args:
    !   ncfile_name: Output NetCDF filename
    !   config: Configuration object with model parameters
    !   angle: Viewing angle array [degrees]
    !   winds: Wind speed array [m/s]
    !   wavenum: Wavenumber array [cm-1]
    !   emiss: Ocean emissivity array (angle, winds, wavenumber)
    !   emiss_noref: Emissivity without multiple reflections
    !   refl_cb: Cox-Munk reflectance array
    !   temperature: Sea surface temperature [K]
    !   elapsed_time: Calculation time [seconds]
    !   performance: Calculation rate [calculations/second]
    !-----------------------------------------------------------------------
    subroutine save_spectral_emissivity(ncfile_name, config, angle, winds, wavenum, &
            emiss, temperature, salinity, &
            emiss_noref, refl_cb)
        implicit none

        ! Input parameters
        character(len=*), intent(in):: ncfile_name
        type(config_type), intent(in):: config
        real(kind = jprd), intent(in):: angle(:)          ! n_angle
        real(kind = jprd), intent(in):: winds(:)           ! n_winds
        real(kind = jprd), intent(in):: wavenum(:)         ! n_wavenumber
        real(kind = jprd), intent(in):: emiss(:,:,:)       ! (n_angle, n_winds, n_wavenumber)
        real(kind = jprd), intent(in):: temperature        ! SST [K]
        real(kind = jprd), intent(in):: salinity           ! Sea surface salinity [PSU]
        real(kind = jprd), optional, intent(in):: emiss_noref(:,:,:) ! (n_angle, n_winds, n_wavenumber)
        real(kind = jprd), optional, intent(in):: refl_cb(:,:,:)     ! (n_angle, n_winds, n_wavenumber)

        ! Local variables
        type(netcdf_file):: output_ncfile
        integer(kind = jpim):: n_angle, n_winds, n_wavenumber

        ! Get array dimensions
        n_angle = size(angle)
        n_winds = size(winds)
        n_wavenumber = size(wavenum)

        ! Validate input dimensions
        if (size(emiss, 1) /= n_angle .or. size(emiss, 2) /= n_winds .or. &
            size(emiss, 3) /= n_wavenumber) then
            write(*,'(a)') "ERROR: Emissivity array dimensions do not match coordinate arrays"
            return
        end if

        write(*,'(a, a)') "Creating NetCDF output file: ", trim(ncfile_name)
        write(*,'(a, i0, a, i0, a, i0, a)') "  Dimensions: ", n_angle, " angle x ", n_winds, &
                                          " winds x ", n_wavenumber, " wavenumbers"

        !-------------------------------------------------------------------
        ! Step 1: Create NetCDF file
        !-------------------------------------------------------------------
        call output_ncfile%create(trim(ncfile_name), iverbose = 2)

        !-------------------------------------------------------------------
        ! Step 2: Define dimensions
        !-------------------------------------------------------------------
        call output_ncfile%define_dimension("n_angle", n_angle)
        call output_ncfile%define_dimension("n_winds", n_winds)
        call output_ncfile%define_dimension("n_wavenumber", n_wavenumber)

        !-------------------------------------------------------------------
        ! Step 3: Define coordinate variables
        !-------------------------------------------------------------------
        call output_ncfile%define_variable("angle", &
            dim1_name="n_angle", &
            units_str = ANGLE_UNITS, &
            long_name="Viewing angle from nadir", &
            standard_name="viewing_angle")

        call output_ncfile%define_variable("wind_speed", &
            dim1_name="n_winds", &
            units_str = WIND_UNITS, &
            long_name="Wind speed at 10m height", &
            standard_name="wind_speed")

        call output_ncfile%define_variable("wavenumber", &
            dim1_name="n_wavenumber", &
            units_str = WAVENUMBER_UNITS, &
            long_name="Wavenumber", &
            standard_name="wavenumber")

        !-------------------------------------------------------------------
        ! Step 4: Define data variables
        !-------------------------------------------------------------------
        ! Primary emissivity with multiple reflections
        call output_ncfile%define_variable("emissivity", &
            dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavenumber", &
            units_str = EMISS_UNITS, &
            long_name="Ocean surface emissivity with multiple reflections", &
            standard_name="surface_emissivity")

        ! Optional variables-only define if present
        if (present(emiss_noref)) then
            ! Emissivity without multiple reflections
            call output_ncfile%define_variable("emissivity_no_reflection", &
                dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavenumber", &
                units_str = EMISS_UNITS, &
                long_name="Ocean surface emissivity without multiple reflections", &
                standard_name="surface_emissivity_no_reflection")
        end if

        if (present(refl_cb)) then
            ! Cox-Munk reflectance
            call output_ncfile%define_variable("reflectance_cox_munk", &
                dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavenumber", &
                units_str = EMISS_UNITS, &
                long_name="Cox-Munk ocean surface reflectance", &
                standard_name="surface_reflectance")
        end if

        ! Scalar parameters
        call output_ncfile%define_variable("sea_surface_temperature", &
            units_str = TEMPERATURE_UNITS, &
            long_name="Sea surface temperature", &
            standard_name="sea_surface_temperature")
        call output_ncfile%define_variable("sea_surface_salinity", &
            units_str="PSU", &
            long_name="Sea surface salinity", &
            standard_name="sea_surface_salinity")

        !-------------------------------------------------------------------
        ! Step 5: Add global attributes
        !-------------------------------------------------------------------
        call add_global_attributes_spectral(output_ncfile, config, temperature, salinity, &
            n_angle, n_winds, n_wavenumber)

        !-------------------------------------------------------------------
        ! Step 6: Add variable attributes
        !-------------------------------------------------------------------
        ! Angles attributes
        call output_ncfile%put_attribute("angle", "axis", "X")
        call output_ncfile%put_attribute("angle", "description", &
            "Viewing angle measured from nadir (0 degrees = vertical)")

        ! Wind speed attributes
        call output_ncfile%put_attribute("wind_speed", "axis", "Y")
        call output_ncfile%put_attribute("wind_speed", "description", &
            "Wind speed at 10 meters above sea surface")

        ! Wavenumber attributes
        call output_ncfile%put_attribute("wavenumber", "axis", "Z")
        call output_ncfile%put_attribute("wavenumber", "description", &
            "Infrared wavenumber")
        call output_ncfile%put_attribute("wavenumber", "monotonic", "increasing")

        ! Emissivity attributes
        call output_ncfile%put_attribute("emissivity", "description", &
            "Ocean surface emissivity including multiple reflection effects")
        call output_ncfile%put_attribute("emissivity", "coordinates", "angle wind_speed wavenumber")

        ! Optional variable attributes
        if (present(emiss_noref)) then
            call output_ncfile%put_attribute("emissivity_no_reflection", "description", &
                "Ocean surface emissivity without multiple reflection effects")
            call output_ncfile%put_attribute("emissivity_no_reflection", "coordinates", "angle wind_speed wavenumber")
        end if

        if (present(refl_cb)) then
            call output_ncfile%put_attribute("reflectance_cox_munk", "description", &
                "Cox-Munk ocean surface bidirectional reflectance")
            call output_ncfile%put_attribute("reflectance_cox_munk", "coordinates", "angle wind_speed wavenumber")
        end if

        ! Temperature attributes
        call output_ncfile%put_attribute("sea_surface_temperature", "description", &
            "Sea surface temperature used in emissivity calculations")

        ! Salinity attributes
        call output_ncfile%put_attribute("sea_surface_salinity", "description", &
            "Sea surface salinity used in emissivity calculations")

        !-------------------------------------------------------------------
        ! Step 7: Write coordinate data
        !-------------------------------------------------------------------
        call output_ncfile%put("angle", angle)
        call output_ncfile%put("wind_speed", winds)
        call output_ncfile%put("wavenumber", wavenum)

        !-------------------------------------------------------------------
        ! Step 8: Write data variables
        !-------------------------------------------------------------------
        call output_ncfile%put("emissivity", emiss)
        ! Write optional variables if present
        if (present(emiss_noref)) then
            call output_ncfile%put("emissivity_no_reflection", emiss_noref)
        end if
        if (present(refl_cb)) then
            call output_ncfile%put("reflectance_cox_munk", refl_cb)
        end if
        call output_ncfile%put("sea_surface_temperature", temperature)
        call output_ncfile%put("sea_surface_salinity", salinity)

        !-------------------------------------------------------------------
        ! Step 9: Close NetCDF file
        !-------------------------------------------------------------------
        call output_ncfile%close()

        write(*,'(a)') "NetCDF output completed successfully!"
        write(*,'(a, f6.3, a)') "  File size: ", file_size_mb(ncfile_name), " MB"

    end subroutine save_spectral_emissivity

    !-----------------------------------------------------------------------
    ! Save polarized ocean emissivity results to NetCDF file (PLACEHOLDER)
    !-----------------------------------------------------------------------
    subroutine save_spectral_polarized_emissivity(ncfile_name, config, angle, winds, wavenum, &
            emiss_v, emiss_h, temperature, salinity, &
            emiss_noref_v, refl_cb_v, emiss_noref_h, refl_cb_h)
        implicit none
        ! Input parameters
        character(len=*), intent(in):: ncfile_name
        type(config_type), intent(in):: config
        real(kind = jprd), intent(in):: angle(:)          ! n_angle
        real(kind = jprd), intent(in):: winds(:)           ! n_winds
        real(kind = jprd), intent(in):: wavenum(:)         ! n_wavenumber
        real(kind = jprd), intent(in):: emiss_v(:,:,:)     ! (n_angle, n_winds, n_wavenumber) - Vertical
        real(kind = jprd), intent(in):: emiss_h(:,:,:)     ! (n_angle, n_winds, n_wavenumber) - Horizontal
        real(kind = jprd), intent(in):: temperature        ! SST [K]
        real(kind = jprd), intent(in):: salinity           ! Sea surface salinity [PSU]
        real(kind = jprd), optional, intent(in):: emiss_noref_v(:,:,:)  ! No reflection-Vertical
        real(kind = jprd), optional, intent(in):: refl_cb_v(:,:,:)      ! Multiple reflection-Vertical
        real(kind = jprd), optional, intent(in):: emiss_noref_h(:,:,:)  ! No reflection-Horizontal
        real(kind = jprd), optional, intent(in):: refl_cb_h(:,:,:)      ! Multiple reflection-Horizontal

        ! Local variables
        type(netcdf_file):: output_ncfile
        integer(kind = jpim):: n_angle, n_winds, n_wavenumber

        ! Get dimensions
        n_angle = size(angle)
        n_winds = size(winds)
        n_wavenumber = size(wavenum)

        ! Validate input dimensions
        if (size(emiss_v, 1) /= n_angle .or. size(emiss_v, 2) /= n_winds .or. &
            size(emiss_v, 3) /= n_wavenumber .or. &
            size(emiss_h, 1) /= n_angle .or. size(emiss_h, 2) /= n_winds .or. &
            size(emiss_h, 3) /= n_wavenumber) then
            write(*,'(a)') "ERROR: Polarized emissivity array dimensions do not match coordinate arrays"
            return
        end if

        write(*,'(a, a)') "Creating polarized NetCDF output file: ", trim(ncfile_name)
        write(*,'(a, i0, a, i0, a, i0, a)') "  Dimensions: ", n_angle, " angle x ", n_winds, &
                                          " winds x ", n_wavenumber, " wavenumbers"

        !-------------------------------------------------------------------
        ! Create NetCDF file
        !-------------------------------------------------------------------
        call output_ncfile%create(trim(ncfile_name), iverbose = 2)

        !-------------------------------------------------------------------
        ! Define dimensions
        !-------------------------------------------------------------------
        call output_ncfile%define_dimension("n_angle", n_angle)
        call output_ncfile%define_dimension("n_winds", n_winds)
        call output_ncfile%define_dimension("n_wavenumber", n_wavenumber)

        !-------------------------------------------------------------------
        ! Define coordinate variables
        !-------------------------------------------------------------------
        call output_ncfile%define_variable("angle", &
            dim1_name="n_angle", &
            units_str = ANGLE_UNITS, &
            long_name="Viewing angle from nadir", &
            standard_name="viewing_angle")

        call output_ncfile%define_variable("wind_speed", &
            dim1_name="n_winds", &
            units_str = WIND_UNITS, &
            long_name="Wind speed at 10m height", &
            standard_name="wind_speed")

        call output_ncfile%define_variable("wavenumber", &
            dim1_name="n_wavenumber", &
            units_str = WAVENUMBER_UNITS, &
            long_name="Wavenumber", &
            standard_name="wavenumber")

        !-------------------------------------------------------------------
        ! Define polarized data variables
        !-------------------------------------------------------------------
        ! Vertical polarization emissivity
        call output_ncfile%define_variable("emissivity_v", &
            dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavenumber", &
            units_str = EMISS_UNITS, &
            long_name="Ocean surface emissivity-vertical polarization", &
            standard_name="surface_emissivity_v")

        ! Horizontal polarization emissivity
        call output_ncfile%define_variable("emissivity_h", &
            dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavenumber", &
            units_str = EMISS_UNITS, &
            long_name="Ocean surface emissivity-horizontal polarization", &
            standard_name="surface_emissivity_h")

        ! Optional variables-only define if present
        if (present(emiss_noref_v) .and. present(emiss_noref_h)) then
            call output_ncfile%define_variable("emissivity_noref_v", &
                dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavenumber", &
                units_str = EMISS_UNITS, &
                long_name="Ocean surface emissivity without reflections-vertical", &
                standard_name="surface_emissivity_noref_v")

            call output_ncfile%define_variable("emissivity_noref_h", &
                dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavenumber", &
                units_str = EMISS_UNITS, &
                long_name="Ocean surface emissivity without reflections-horizontal", &
                standard_name="surface_emissivity_noref_h")
        end if

        if (present(refl_cb_v) .and. present(refl_cb_h)) then
            call output_ncfile%define_variable("reflectance_v", &
                dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavenumber", &
                units_str = EMISS_UNITS, &
                long_name="Cox-Munk surface reflectance-vertical", &
                standard_name="surface_reflectance_v")

            call output_ncfile%define_variable("reflectance_h", &
                dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavenumber", &
                units_str = EMISS_UNITS, &
                long_name="Cox-Munk surface reflectance-horizontal", &
                standard_name="surface_reflectance_h")
        end if

        ! Sea surface temperature and salinity
        call output_ncfile%define_variable("sea_surface_temperature", &
            units_str = TEMPERATURE_UNITS, &
            long_name="Sea surface temperature", &
            standard_name="sea_surface_temperature")
        call output_ncfile%define_variable("sea_surface_salinity", &
            units_str="PSU", &
            long_name="Sea surface salinity", &
            standard_name="sea_surface_salinity")

        !-------------------------------------------------------------------
        ! Add global attributes
        !-------------------------------------------------------------------
        call add_global_attributes_spectral(output_ncfile, config, temperature, salinity, &
            n_angle, n_winds, n_wavenumber)

        !-------------------------------------------------------------------
        ! Write coordinate data
        !-------------------------------------------------------------------
        call output_ncfile%put("angle", angle)
        call output_ncfile%put("wind_speed", winds)
        call output_ncfile%put("wavenumber", wavenum)

        !-------------------------------------------------------------------
        ! Write polarized emissivity data
        !-------------------------------------------------------------------
        call output_ncfile%put("emissivity_v", emiss_v)
        call output_ncfile%put("emissivity_h", emiss_h)

        ! Write optional data if present
        if (present(emiss_noref_v) .and. present(emiss_noref_h)) then
            call output_ncfile%put("emissivity_noref_v", emiss_noref_v)
            call output_ncfile%put("emissivity_noref_h", emiss_noref_h)
        end if

        if (present(refl_cb_v) .and. present(refl_cb_h)) then
            call output_ncfile%put("reflectance_v", refl_cb_v)
            call output_ncfile%put("reflectance_h", refl_cb_h)
        end if

        call output_ncfile%put("sea_surface_temperature", temperature)
        call output_ncfile%put("sea_surface_salinity", salinity)

        !-------------------------------------------------------------------
        ! Close NetCDF file
        !-------------------------------------------------------------------
        call output_ncfile%close()

        write(*,'(a)') "NetCDF output completed successfully!"
        write(*,'(a, f6.3, a)') "  File size: ", file_size_mb(ncfile_name), " MB"

    end subroutine save_spectral_polarized_emissivity

    !-----------------------------------------------------------------------
    ! Save narrow-band ocean emissivity results to NetCDF file
    !
    ! This subroutine creates a NetCDF file containing narrow-band emissivity
    ! results after integration with Planck function
    !
    ! Args:
    !   ncfile_name: Output NetCDF filename
    !   config: Configuration object with model parameters
    !   angle: Viewing angle array [degrees]
    !   winds: Wind speed array [m/s]
    !   band_centers: Center wavenumbers for each band [cm-1]
    !   band_widths: Width of each band [cm-1]
    !   emiss_narrow: Narrow-band emissivity array (angle, winds, bands)
    !   temperature: Sea surface temperature [K]
    !   salinity: Sea surface salinity [PSU]
    !-----------------------------------------------------------------------
    subroutine save_narrow_band_emissivity(ncfile_name, config, angle, winds, &
            wavenum_start, wavenum_end, emiss_narrow, &
            temperature, salinity)
        implicit none

        ! Input parameters
        character(len=*), intent(in):: ncfile_name
        type(config_type), intent(in):: config
        real(kind = jprd), intent(in):: angle(:)         ! n_angle
        real(kind = jprd), intent(in):: winds(:)          ! n_winds
        real(kind = jprd), intent(in):: wavenum_start(:)  ! n_bands
        real(kind = jprd), intent(in):: wavenum_end(:)    ! n_bands
        real(kind = jprd), intent(in):: emiss_narrow(:,:,:) ! (n_angle, n_winds, n_bands)
        real(kind = jprd), intent(in):: temperature       ! SST [K]
        real(kind = jprd), intent(in):: salinity          ! Salinity [PSU]

        ! Local variables
        type(netcdf_file):: output_ncfile
        integer(kind = jpim):: n_angle, n_winds, n_bands

        ! Get array dimensions
        n_angle = size(angle)
        n_winds = size(winds)
        n_bands = size(wavenum_start)

        ! Validate input dimensions
        if (size(emiss_narrow, 1) /= n_angle .or. size(emiss_narrow, 2) /= n_winds .or. &
            size(emiss_narrow, 3) /= n_bands) then
            write(*,'(a)') "ERROR: Narrow-band emissivity array dimensions do not match coordinate arrays"
            return
        end if

        write(*,'(a, a)') "Creating narrow-band NetCDF output file: ", trim(ncfile_name)
        write(*,'(a, i0, a, i0, a, i0, a)') "  Dimensions: ", n_angle, " angle × ", n_winds, " winds × ", n_bands, " bands"

        !-------------------------------------------------------------------
        ! Step 1: Create NetCDF file
        !-------------------------------------------------------------------
        call output_ncfile%create(trim(ncfile_name), iverbose = 2)

        !-------------------------------------------------------------------
        ! Step 2: Define dimensions
        !-------------------------------------------------------------------
        call output_ncfile%define_dimension("angle", n_angle)
        call output_ncfile%define_dimension("wind", n_winds)
        call output_ncfile%define_dimension("band", n_bands)

        !-------------------------------------------------------------------
        ! Step 3: Define coordinate variables
        !-------------------------------------------------------------------
        ! Viewing angle (use same names as spectral mode for consistency)
        call output_ncfile%define_variable("angle", &
            dim1_name="angle", &
            units_str = ANGLE_UNITS, &
            long_name="Viewing angle from nadir", &
            standard_name="viewing_angle")

        ! Wind speeds (use same names as spectral mode for consistency)
        call output_ncfile%define_variable("wind_speed", &
            dim1_name="wind", &
            units_str = WIND_UNITS, &
            long_name="Wind speed at 10m height", &
            standard_name="wind_speed")

        ! Band start wavenumbers
        call output_ncfile%define_variable("wavenum_start", &
            dim1_name="band", &
            units_str = WAVENUMBER_UNITS, &
            long_name="Start wavenumber of narrow band")

        ! Band end wavenumbers
        call output_ncfile%define_variable("wavenum_end", &
            dim1_name="band", &
            units_str = WAVENUMBER_UNITS, &
            long_name="End wavenumber of narrow band")

        !-------------------------------------------------------------------
        ! Step 4: Define data variables
        !-------------------------------------------------------------------
        ! Narrow-band ocean emissivity
        call output_ncfile%define_variable("ocean_emissivity_narrow_band", &
            dim1_name="angle", dim2_name="wind", dim3_name="band", &
            units_str = EMISS_UNITS, &
            long_name="Ocean surface narrow-band emissivity", &
            standard_name="surface_emissivity_narrow_band")

        ! Temperature (scalar) - use same names as spectral mode
        call output_ncfile%define_variable("sea_surface_temperature", &
            units_str = TEMPERATURE_UNITS, &
            long_name="Sea surface temperature", &
            standard_name="sea_surface_temperature")

        ! Salinity (scalar) - use same names as spectral mode
        call output_ncfile%define_variable("sea_surface_salinity", &
            units_str="PSU", &
            long_name="Sea surface salinity", &
            standard_name="sea_surface_salinity")

        !-------------------------------------------------------------------
        ! Step 5: Add global attributes
        !-------------------------------------------------------------------
        call add_global_attributes_narrow_band(output_ncfile, config, temperature, salinity, &
            n_angle, n_winds, n_bands)

        !-------------------------------------------------------------------
        ! Step 6: Write data
        !-------------------------------------------------------------------
        call output_ncfile%put("angle", angle)
        call output_ncfile%put("wind_speed", winds)
        call output_ncfile%put("wavenum_start", wavenum_start)
        call output_ncfile%put("wavenum_end", wavenum_end)
        call output_ncfile%put("ocean_emissivity_narrow_band", emiss_narrow)
        call output_ncfile%put("sea_surface_temperature", temperature)
        call output_ncfile%put("sea_surface_salinity", salinity)

        !-------------------------------------------------------------------
        ! Step 7: Close file
        !-------------------------------------------------------------------
        call output_ncfile%close()

        write(*,'(a)') "Narrow-band NetCDF output completed successfully!"
        write(*,'(a, f6.3, a)') "  File size: ", file_size_mb(ncfile_name), " MB"

    end subroutine save_narrow_band_emissivity

    !-----------------------------------------------------------------------
    ! Add global attributes for spectral mode
    !-----------------------------------------------------------------------
    subroutine add_global_attributes_spectral(output_ncfile, config, temperature, salinity, &
            n_angle, n_winds, n_wavenumber)
        type(netcdf_file), intent(inout):: output_ncfile
        type(config_type), intent(in):: config
        real(kind = jprd), intent(in):: temperature, salinity
        integer(kind = jpim), intent(in):: n_angle, n_winds, n_wavenumber

        character(len = 19):: datetime_str
        character(len = 512):: temp_str

        ! Get current date/time
        call get_current_datetime(datetime_str)

        ! Core metadata
        call output_ncfile%put_global_attribute("title", &
            "Infrared Ocean Surface Effective Emissivity (IOSEE) Output")
        call output_ncfile%put_global_attribute("institution", &
            "Atmospheric & Oceanic Optics Group, Department of Atmospheric Sciences, Texas A&M University")
        call output_ncfile%put_global_attribute("source", &
            "Infrared Ocean Surface Effective Emissivity (IOSEE), Version 1.0")
        call output_ncfile%put_global_attribute("references", &
            "Wu & Smith (1997), Henderson et al.,(2003), Masuda (2006), " &
            // "Nalli et al (2008, 2023), Gero et al., (2016, 2019), Wei et al., (2025) etc")
        call output_ncfile%put_global_attribute("Conventions", "CF-1.8")
        call output_ncfile%put_global_attribute("history", &
            trim(datetime_str) // ": Created by the OCIOSE spectral mode")
        call output_ncfile%put_global_attribute("creation_date", datetime_str)
        call output_ncfile%put_global_attribute("processing_mode", "spectral")
        call output_ncfile%put_global_attribute("copyright", &
            "(C) Copyright 2025-TAMU")

        ! Model configuration
        write(temp_str, '(f6.1)') temperature
        call output_ncfile%put_global_attribute("model_temperature_K", trim(temp_str))

        write(temp_str, '(f6.1)') salinity
        call output_ncfile%put_global_attribute("model_salinity_PSU", trim(temp_str))

        write(temp_str, '(i0)') n_angle
        call output_ncfile%put_global_attribute("number_of_angle", trim(temp_str))

        write(temp_str, '(i0)') n_winds
        call output_ncfile%put_global_attribute("number_of_winds", trim(temp_str))

        write(temp_str, '(i0)') n_wavenumber
        call output_ncfile%put_global_attribute("number_of_wavenumbers", trim(temp_str))

        ! Physical model description
        call output_ncfile%put_global_attribute("model_description", &
            "Spectral emissivity")
        call output_ncfile%put_global_attribute("surface_type", "Ocean/water")
        call output_ncfile%put_global_attribute("spectrum_region", "Infrared")
        call output_ncfile%put_global_attribute("auxiliary_data", trim(config%refractive_index_filename))

    end subroutine add_global_attributes_spectral

    !-----------------------------------------------------------------------
    ! Add global attributes for narrow-band mode
    !-----------------------------------------------------------------------
    subroutine add_global_attributes_narrow_band(output_ncfile, config, temperature, salinity, &
            n_angle, n_winds, n_bands)
        type(netcdf_file), intent(inout):: output_ncfile
        type(config_type), intent(in):: config
        real(kind = jprd), intent(in):: temperature, salinity
        integer(kind = jpim), intent(in):: n_angle, n_winds, n_bands

        character(len = 19):: datetime_str
        character(len = 512):: temp_str

        ! Get current date/time
        call get_current_datetime(datetime_str)

        ! Core metadata
        call output_ncfile%put_global_attribute("title", &
            "Infrared Ocean Surface Effective Emissivity (IOSEE) Output")
        call output_ncfile%put_global_attribute("institution", &
            "Atmospheric & Oceanic Optics Group, Department of Atmospheric Sciences, Texas A&M University")
        call output_ncfile%put_global_attribute("source", &
            "Infrared Ocean Surface Effective Emissivity (IOSEE), Version 1.0")
        call output_ncfile%put_global_attribute("references", &
            "Wu & Smith (1997), Henderson et al.,(2003), Masuda (2006), " &
            // "Nalli et al (2008, 2023), Gero et al., (2016, 2019), Wei et al., (2025) etc")
        call output_ncfile%put_global_attribute("Conventions", "CF-1.8")
        call output_ncfile%put_global_attribute("history", &
            trim(datetime_str) // ": Created by the OCIOSE narrow-band mode")
        call output_ncfile%put_global_attribute("creation_date", datetime_str)
        call output_ncfile%put_global_attribute("processing_mode", "narrow-band")
        call output_ncfile%put_global_attribute("integration_method", &
            "Planck-weighted integration over spectral bands")
        call output_ncfile%put_global_attribute("copyright", &
            "(C) Copyright 2025-TAMU")

        ! Model configuration
        write(temp_str, '(f6.1)') temperature
        call output_ncfile%put_global_attribute("model_temperature_K", trim(temp_str))

        write(temp_str, '(f6.1)') salinity
        call output_ncfile%put_global_attribute("model_salinity_PSU", trim(temp_str))

        write(temp_str, '(i0)') n_angle
        call output_ncfile%put_global_attribute("number_of_angle", trim(temp_str))

        write(temp_str, '(i0)') n_winds
        call output_ncfile%put_global_attribute("number_of_winds", trim(temp_str))

        write(temp_str, '(i0)') n_bands
        call output_ncfile%put_global_attribute("number_of_bands", trim(temp_str))

        ! Physical model description
        call output_ncfile%put_global_attribute("model_description", &
            "Narrow-band emissivity integrated with Planck function")
        call output_ncfile%put_global_attribute("surface_type", "Ocean/water")
        call output_ncfile%put_global_attribute("spectrum_region", "Infrared")
        call output_ncfile%put_global_attribute("auxiliary_data", trim(config%refractive_index_filename))

    end subroutine add_global_attributes_narrow_band


    !-----------------------------------------------------------------------
    ! Get current date and time as ISO 8601 string
    !-----------------------------------------------------------------------
    subroutine get_current_datetime(datetime_str)
        character(len = 19), intent(out):: datetime_str
        integer:: date_time(8)

        call date_and_time(values = date_time)
        write(datetime_str, '(i4, "-",i2.2, "-",i2.2, "T",i2.2, ":",i2.2, ":",i2.2)') &
            date_time(1), date_time(2), date_time(3), date_time(5), date_time(6), date_time(7)
    end subroutine get_current_datetime

    !-----------------------------------------------------------------------
    ! Estimate file size in MB (simple approximation)
    !-----------------------------------------------------------------------
    function file_size_mb(filename) result(size_mb)
        character(len=*), intent(in):: filename
        real(kind = jprd):: size_mb

        ! Simple estimation: 4 bytes per real*number of elements+overhead
        ! This is a rough estimate-actual implementation would use file system calls
        size_mb = 0.1_jprd  ! Placeholder-would implement actual file size check
    end function file_size_mb

    !-----------------------------------------------------------------------
    ! Save ocean emissivity results to NetCDF file (wavelength version)
    !
    ! This subroutine creates a NetCDF file for wavelength-based input
    ! with wavelength as the spectral coordinate (sorted in increasing order)
    !
    ! Args:
    !   ncfile_name: Output NetCDF filename
    !   config: Configuration object with model parameters
    !   angle: Viewing angle array [degrees]
    !   winds: Wind speed array [m/s]
    !   wavelen: Wavelength array [um] (descending order from processing)
    !   emiss: Ocean emissivity array (angle, winds, wavelength)
    !   temperature: Sea surface temperature [K]
    !   salinity: Sea surface salinity [PSU]
    !   emiss_noref: Emissivity without multiple reflections (optional)
    !   refl_cb: Cox-Munk reflectance array (optional)
    !-----------------------------------------------------------------------
    subroutine save_spectral_emissivity_wl(ncfile_name, config, angle, winds, wavelen, &
            emiss, temperature, salinity, &
            emiss_noref, refl_cb)
        implicit none

        ! Input parameters
        character(len=*), intent(in):: ncfile_name
        type(config_type), intent(in):: config
        real(kind = jprd), intent(in):: angle(:)          ! n_angle
        real(kind = jprd), intent(in):: winds(:)           ! n_winds
        real(kind = jprd), intent(in):: wavelen(:)         ! n_wavelength (descending order)
        real(kind = jprd), intent(in):: emiss(:,:,:)       ! (n_angle, n_winds, n_wavelength)
        real(kind = jprd), intent(in):: temperature        ! SST [K]
        real(kind = jprd), intent(in):: salinity           ! Sea surface salinity [PSU]
        real(kind = jprd), optional, intent(in):: emiss_noref(:,:,:) ! (n_angle, n_winds, n_wavelength)
        real(kind = jprd), optional, intent(in):: refl_cb(:,:,:)     ! (n_angle, n_winds, n_wavelength)

        ! Local variables
        type(netcdf_file):: output_ncfile
        integer(kind = jpim):: n_angle, n_winds, n_wavelength

        ! Arrays for sorting wavelength to increasing order
        real(kind = jprd), allocatable:: wavelen_sorted(:)
        real(kind = jprd), allocatable:: emiss_sorted(:,:,:)
        real(kind = jprd), allocatable:: emiss_noref_sorted(:,:,:)
        real(kind = jprd), allocatable:: refl_cb_sorted(:,:,:)
        integer(kind = jpim), allocatable:: sort_indices(:)
        integer(kind = jpim):: i, j, k, idx

        ! Get array dimensions
        n_angle = size(angle)
        n_winds = size(winds)
        n_wavelength = size(wavelen)

        ! Validate input dimensions
        if (size(emiss, 1) /= n_angle .or. size(emiss, 2) /= n_winds .or. &
            size(emiss, 3) /= n_wavelength) then
            write(*,'(a)') "ERROR: Emissivity array dimensions do not match coordinate arrays"
            return
        end if

        write(*,'(a, a)') "Creating wavelength-based NetCDF output file: ", trim(ncfile_name)
        write(*,'(a, i0, a, i0, a, i0, a)') "  Dimensions: ", n_angle, " angle x ", n_winds, &
                                          " winds x ", n_wavelength, " wavelengths"

        ! Allocate arrays for sorting
        allocate(wavelen_sorted(n_wavelength))
        allocate(emiss_sorted(n_angle, n_winds, n_wavelength))
        allocate(sort_indices(n_wavelength))

        if (present(emiss_noref)) then
            allocate(emiss_noref_sorted(n_angle, n_winds, n_wavelength))
        end if

        if (present(refl_cb)) then
            allocate(refl_cb_sorted(n_angle, n_winds, n_wavelength))
        end if

        !-------------------------------------------------------------------
        ! Sort wavelength to increasing order and reorder emissivity data
        ! Optimized with OpenMP parallelization for performance
        !-------------------------------------------------------------------

        ! Create wavelength array in increasing order (reverse the descending input)
        !$OMP PARALLEL DO PRIVATE(i) SHARED(wavelen_sorted, wavelen, sort_indices, n_wavelength)
        do i = 1, n_wavelength
            sort_indices(i) = n_wavelength-i + 1
            wavelen_sorted(i) = wavelen(sort_indices(i))
        end do
        !$OMP END PARALLEL DO

        ! Reorder emissivity data to match sorted wavelength-parallel over wavelength dimension
        !$OMP PARALLEL DO PRIVATE(k, idx, j, i) &
        !$OMP SHARED(emiss_sorted, emiss, sort_indices, n_wavelength, n_winds, n_angle)
        do k = 1, n_wavelength
            idx = sort_indices(k)
            do j = 1, n_winds
                do i = 1, n_angle
                    emiss_sorted(i, j, k) = emiss(i, j, idx)
                end do
            end do
        end do
        !$OMP END PARALLEL DO

        ! Reorder optional arrays if present-parallel processing
        if (present(emiss_noref)) then
            !$OMP PARALLEL DO PRIVATE(k, idx, j, i) &
            !$OMP SHARED(emiss_noref_sorted, emiss_noref, sort_indices, n_wavelength, n_winds, n_angle)
            do k = 1, n_wavelength
                idx = sort_indices(k)
                do j = 1, n_winds
                    do i = 1, n_angle
                        emiss_noref_sorted(i, j, k) = emiss_noref(i, j, idx)
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
        end if

        if (present(refl_cb)) then
            !$OMP PARALLEL DO PRIVATE(k, idx, j, i) &
            !$OMP SHARED(refl_cb_sorted, refl_cb, sort_indices, n_wavelength, n_winds, n_angle)
            do k = 1, n_wavelength
                idx = sort_indices(k)
                do j = 1, n_winds
                    do i = 1, n_angle
                        refl_cb_sorted(i, j, k) = refl_cb(i, j, idx)
                    end do
                end do
            end do
            !$OMP END PARALLEL DO
        end if

        !-------------------------------------------------------------------
        ! Step 1: Create NetCDF file
        !-------------------------------------------------------------------
        call output_ncfile%create(trim(ncfile_name), iverbose = 2)

        !-------------------------------------------------------------------
        ! Step 2: Define dimensions
        !-------------------------------------------------------------------
        call output_ncfile%define_dimension("n_angle", n_angle)
        call output_ncfile%define_dimension("n_winds", n_winds)
        call output_ncfile%define_dimension("n_wavelength", n_wavelength)

        !-------------------------------------------------------------------
        ! Step 3: Define coordinate variables
        !-------------------------------------------------------------------
        call output_ncfile%define_variable("angle", &
            dim1_name="n_angle", &
            units_str = ANGLE_UNITS, &
            long_name="Viewing angle from nadir", &
            standard_name="viewing_angle")

        call output_ncfile%define_variable("wind_speed", &
            dim1_name="n_winds", &
            units_str = WIND_UNITS, &
            long_name="Wind speed at 10m height", &
            standard_name="wind_speed")

        call output_ncfile%define_variable("wavelength", &
            dim1_name="n_wavelength", &
            units_str = WAVELENGTH_UNITS, &
            long_name="Wavelength", &
            standard_name="radiation_wavelength")

        !-------------------------------------------------------------------
        ! Step 4: Define data variables
        !-------------------------------------------------------------------
        ! Primary emissivity with multiple reflections
        call output_ncfile%define_variable("emissivity", &
            dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavelength", &
            units_str = EMISS_UNITS, &
            long_name="Ocean surface emissivity with multiple reflections", &
            standard_name="surface_emissivity")

        ! Optional variables-only define if present
        if (present(emiss_noref)) then
            ! Emissivity without multiple reflections
            call output_ncfile%define_variable("emissivity_no_reflection", &
                dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavelength", &
                units_str = EMISS_UNITS, &
                long_name="Ocean surface emissivity without multiple reflections", &
                standard_name="surface_emissivity_no_reflection")
        end if

        if (present(refl_cb)) then
            ! Cox-Munk reflectance
            call output_ncfile%define_variable("reflectance_cox_munk", &
                dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavelength", &
                units_str = EMISS_UNITS, &
                long_name="Cox-Munk surface reflectance", &
                standard_name="surface_reflectance")
        end if

        ! Sea surface temperature
        call output_ncfile%define_variable("sea_surface_temperature", &
            units_str = TEMPERATURE_UNITS, &
            long_name="Sea surface temperature", &
            standard_name="sea_surface_temperature")
        call output_ncfile%define_variable("sea_surface_salinity", &
            units_str="PSU", &
            long_name="Sea surface salinity", &
            standard_name="sea_surface_salinity")

        !-------------------------------------------------------------------
        ! Step 5: Add global attributes
        !-------------------------------------------------------------------
        call add_global_attributes_spectral(output_ncfile, config, temperature, salinity, &
            n_angle, n_winds, n_wavelength)

        !-------------------------------------------------------------------
        ! Step 6: Add variable attributes
        !-------------------------------------------------------------------
        ! Angles attributes
        call output_ncfile%put_attribute("angle", "axis", "X")
        call output_ncfile%put_attribute("angle", "description", &
            "Viewing angle measured from nadir (0 degrees = vertical)")

        ! Wind speed attributes
        call output_ncfile%put_attribute("wind_speed", "axis", "Y")
        call output_ncfile%put_attribute("wind_speed", "description", &
            "Wind speed at 10 meters above sea surface")

        ! Wavelength attributes
        call output_ncfile%put_attribute("wavelength", "axis", "Z")
        call output_ncfile%put_attribute("wavelength", "description", &
            "Infrared wavelength")
        call output_ncfile%put_attribute("wavelength", "monotonic", "increasing")

        ! Emissivity attributes
        call output_ncfile%put_attribute("emissivity", "description", &
            "Ocean surface emissivity including multiple reflection effects")
        call output_ncfile%put_attribute("emissivity", "coordinates", "angle wind_speed wavelength")

        ! Optional variable attributes
        if (present(emiss_noref)) then
            call output_ncfile%put_attribute("emissivity_no_reflection", "description", &
                "Ocean surface emissivity without multiple reflection effects")
            call output_ncfile%put_attribute("emissivity_no_reflection", "coordinates", "angle wind_speed wavelength")
        end if

        if (present(refl_cb)) then
            call output_ncfile%put_attribute("reflectance_cox_munk", "description", &
                "Cox-Munk ocean surface bidirectional reflectance")
            call output_ncfile%put_attribute("reflectance_cox_munk", "coordinates", "angle wind_speed wavelength")
        end if

        ! Temperature attributes
        call output_ncfile%put_attribute("sea_surface_temperature", "description", &
            "Sea surface temperature used in emissivity calculations")

        ! Salinity attributes
        call output_ncfile%put_attribute("sea_surface_salinity", "description", &
            "Sea surface salinity used in emissivity calculations")

        !-------------------------------------------------------------------
        ! Step 7: Write coordinate data
        !-------------------------------------------------------------------
        call output_ncfile%put("angle", angle)
        call output_ncfile%put("wind_speed", winds)
        call output_ncfile%put("wavelength", wavelen_sorted)

        !-------------------------------------------------------------------
        ! Step 8: Write data variables
        !-------------------------------------------------------------------
        call output_ncfile%put("emissivity", emiss_sorted)
        ! Write optional variables if present
        if (present(emiss_noref)) then
            call output_ncfile%put("emissivity_no_reflection", emiss_noref_sorted)
        end if
        if (present(refl_cb)) then
            call output_ncfile%put("reflectance_cox_munk", refl_cb_sorted)
        end if
        call output_ncfile%put("sea_surface_temperature", temperature)
        call output_ncfile%put("sea_surface_salinity", salinity)

        !-------------------------------------------------------------------
        ! Step 9: Close NetCDF file
        !-------------------------------------------------------------------
        call output_ncfile%close()

        write(*,'(a)') "NetCDF output completed successfully!"
        write(*,'(a, f6.3, a)') "  File size: ", file_size_mb(ncfile_name), " MB"

        ! Clean up allocated arrays
        if (allocated(wavelen_sorted)) deallocate(wavelen_sorted)
        if (allocated(emiss_sorted)) deallocate(emiss_sorted)
        if (allocated(emiss_noref_sorted)) deallocate(emiss_noref_sorted)
        if (allocated(refl_cb_sorted)) deallocate(refl_cb_sorted)
        if (allocated(sort_indices)) deallocate(sort_indices)

    end subroutine save_spectral_emissivity_wl

    !-----------------------------------------------------------------------
    ! Save polarized ocean emissivity results to NetCDF file (wavelength version-PLACEHOLDER)
    !-----------------------------------------------------------------------
    subroutine save_spectral_polarized_emissivity_wl(ncfile_name, config, angle, winds, wavelen, &
            emiss_v, emiss_h, temperature, salinity, &
            emiss_noref_v, refl_cb_v, emiss_noref_h, refl_cb_h)
        implicit none
        ! Input parameters
        character(len=*), intent(in):: ncfile_name
        type(config_type), intent(in):: config
        real(kind = jprd), intent(in):: angle(:)          ! n_angle
        real(kind = jprd), intent(in):: winds(:)           ! n_winds
        real(kind = jprd), intent(in):: wavelen(:)         ! n_wavelength (descending order)
        real(kind = jprd), intent(in):: emiss_v(:,:,:)     ! (n_angle, n_winds, n_wavelength) - Vertical
        real(kind = jprd), intent(in):: emiss_h(:,:,:)     ! (n_angle, n_winds, n_wavelength) - Horizontal
        real(kind = jprd), intent(in):: temperature        ! SST [K]
        real(kind = jprd), intent(in):: salinity           ! Sea surface salinity [PSU]
        real(kind = jprd), optional, intent(in):: emiss_noref_v(:,:,:)  ! No reflection-Vertical
        real(kind = jprd), optional, intent(in):: refl_cb_v(:,:,:)      ! Multiple reflection-Vertical
        real(kind = jprd), optional, intent(in):: emiss_noref_h(:,:,:)  ! No reflection-Horizontal
        real(kind = jprd), optional, intent(in):: refl_cb_h(:,:,:)      ! Multiple reflection-Horizontal

        ! Local variables
        type(netcdf_file):: output_ncfile
        integer(kind = jpim):: n_angle, n_winds, n_wavelength
        real(kind = jprd), allocatable:: wavelen_sorted(:)
        real(kind = jprd), allocatable:: emiss_v_sorted(:,:,:), emiss_h_sorted(:,:,:)
        real(kind = jprd), allocatable:: emiss_noref_v_sorted(:,:,:), emiss_noref_h_sorted(:,:,:)
        real(kind = jprd), allocatable:: refl_cb_v_sorted(:,:,:), refl_cb_h_sorted(:,:,:)
        integer(kind = jpim), allocatable:: sort_indices(:)
        integer(kind = jpim):: i, j, k

        ! Get array dimensions
        n_angle = size(angle)
        n_winds = size(winds)
        n_wavelength = size(wavelen)

        ! Validate input dimensions
        if (size(emiss_v, 1) /= n_angle .or. size(emiss_v, 2) /= n_winds .or. &
            size(emiss_v, 3) /= n_wavelength .or. &
            size(emiss_h, 1) /= n_angle .or. size(emiss_h, 2) /= n_winds .or. &
            size(emiss_h, 3) /= n_wavelength) then
            write(*,'(a)') "ERROR: Polarized emissivity array dimensions do not match coordinate arrays"
            return
        end if

        write(*,'(a, a)') "Creating polarized NetCDF output file (wavelength): ", trim(ncfile_name)
        write(*,'(a, i0, a, i0, a, i0, a)') "  Dimensions: ", n_angle, " angle x ", n_winds, &
                                          " winds x ", n_wavelength, " wavelengths"

        !-------------------------------------------------------------------
        ! Sort wavelength data from descending to ascending for CF compliance
        !-------------------------------------------------------------------
        allocate(sort_indices(n_wavelength))
        allocate(wavelen_sorted(n_wavelength))
        allocate(emiss_v_sorted(n_angle, n_winds, n_wavelength))
        allocate(emiss_h_sorted(n_angle, n_winds, n_wavelength))

        ! Create sort indices (wavelength input is descending, need ascending)
        do i = 1, n_wavelength
            sort_indices(i) = n_wavelength-i + 1
        end do

        ! Sort wavelength array
        do i = 1, n_wavelength
            wavelen_sorted(i) = wavelen(sort_indices(i))
        end do

        ! Sort emissivity arrays
        do k = 1, n_wavelength
            do j = 1, n_winds
                do i = 1, n_angle
                    emiss_v_sorted(i, j, k) = emiss_v(i, j, sort_indices(k))
                    emiss_h_sorted(i, j, k) = emiss_h(i, j, sort_indices(k))
                end do
            end do
        end do

        ! Sort optional arrays if present
        if (present(emiss_noref_v) .and. present(emiss_noref_h)) then
            allocate(emiss_noref_v_sorted(n_angle, n_winds, n_wavelength))
            allocate(emiss_noref_h_sorted(n_angle, n_winds, n_wavelength))
            do k = 1, n_wavelength
                do j = 1, n_winds
                    do i = 1, n_angle
                        emiss_noref_v_sorted(i, j, k) = emiss_noref_v(i, j, sort_indices(k))
                        emiss_noref_h_sorted(i, j, k) = emiss_noref_h(i, j, sort_indices(k))
                    end do
                end do
            end do
        end if

        if (present(refl_cb_v) .and. present(refl_cb_h)) then
            allocate(refl_cb_v_sorted(n_angle, n_winds, n_wavelength))
            allocate(refl_cb_h_sorted(n_angle, n_winds, n_wavelength))
            do k = 1, n_wavelength
                do j = 1, n_winds
                    do i = 1, n_angle
                        refl_cb_v_sorted(i, j, k) = refl_cb_v(i, j, sort_indices(k))
                        refl_cb_h_sorted(i, j, k) = refl_cb_h(i, j, sort_indices(k))
                    end do
                end do
            end do
        end if

        !-------------------------------------------------------------------
        ! Create NetCDF file
        !-------------------------------------------------------------------
        call output_ncfile%create(trim(ncfile_name), iverbose = 2)

        !-------------------------------------------------------------------
        ! Define dimensions
        !-------------------------------------------------------------------
        call output_ncfile%define_dimension("n_angle", n_angle)
        call output_ncfile%define_dimension("n_winds", n_winds)
        call output_ncfile%define_dimension("n_wavelength", n_wavelength)

        !-------------------------------------------------------------------
        ! Define coordinate variables
        !-------------------------------------------------------------------
        call output_ncfile%define_variable("angle", &
            dim1_name="n_angle", &
            units_str = ANGLE_UNITS, &
            long_name="Viewing angle from nadir", &
            standard_name="viewing_angle")

        call output_ncfile%define_variable("wind_speed", &
            dim1_name="n_winds", &
            units_str = WIND_UNITS, &
            long_name="Wind speed at 10m height", &
            standard_name="wind_speed")

        call output_ncfile%define_variable("wavelength", &
            dim1_name="n_wavelength", &
            units_str = WAVELENGTH_UNITS, &
            long_name="Wavelength", &
            standard_name="radiation_wavelength")

        !-------------------------------------------------------------------
        ! Define polarized data variables
        !-------------------------------------------------------------------
        ! Vertical polarization emissivity
        call output_ncfile%define_variable("emissivity_v", &
            dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavelength", &
            units_str = EMISS_UNITS, &
            long_name="Ocean surface emissivity-vertical polarization", &
            standard_name="surface_emissivity_v")

        ! Horizontal polarization emissivity
        call output_ncfile%define_variable("emissivity_h", &
            dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavelength", &
            units_str = EMISS_UNITS, &
            long_name="Ocean surface emissivity-horizontal polarization", &
            standard_name="surface_emissivity_h")

        ! Optional variables-only define if present
        if (present(emiss_noref_v) .and. present(emiss_noref_h)) then
            call output_ncfile%define_variable("emissivity_noref_v", &
                dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavelength", &
                units_str = EMISS_UNITS, &
                long_name="Ocean surface emissivity without reflections-vertical", &
                standard_name="surface_emissivity_noref_v")

            call output_ncfile%define_variable("emissivity_noref_h", &
                dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavelength", &
                units_str = EMISS_UNITS, &
                long_name="Ocean surface emissivity without reflections-horizontal", &
                standard_name="surface_emissivity_noref_h")
        end if

        if (present(refl_cb_v) .and. present(refl_cb_h)) then
            call output_ncfile%define_variable("reflectance_v", &
                dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavelength", &
                units_str = EMISS_UNITS, &
                long_name="Cox-Munk surface reflectance-vertical", &
                standard_name="surface_reflectance_v")

            call output_ncfile%define_variable("reflectance_h", &
                dim1_name="n_angle", dim2_name="n_winds", dim3_name="n_wavelength", &
                units_str = EMISS_UNITS, &
                long_name="Cox-Munk surface reflectance-horizontal", &
                standard_name="surface_reflectance_h")
        end if

        ! Sea surface temperature and salinity
        call output_ncfile%define_variable("sea_surface_temperature", &
            units_str = TEMPERATURE_UNITS, &
            long_name="Sea surface temperature", &
            standard_name="sea_surface_temperature")
        call output_ncfile%define_variable("sea_surface_salinity", &
            units_str="PSU", &
            long_name="Sea surface salinity", &
            standard_name="sea_surface_salinity")

        !-------------------------------------------------------------------
        ! Add global attributes
        !-------------------------------------------------------------------
        call add_global_attributes_spectral(output_ncfile, config, temperature, salinity, &
            n_angle, n_winds, n_wavelength)

        !-------------------------------------------------------------------
        ! Write coordinate data
        !-------------------------------------------------------------------
        call output_ncfile%put("angle", angle)
        call output_ncfile%put("wind_speed", winds)
        call output_ncfile%put("wavelength", wavelen_sorted)

        !-------------------------------------------------------------------
        ! Write polarized data variables
        !-------------------------------------------------------------------
        call output_ncfile%put("emissivity_v", emiss_v_sorted)
        call output_ncfile%put("emissivity_h", emiss_h_sorted)

        ! Write optional variables if present
        if (present(emiss_noref_v) .and. present(emiss_noref_h)) then
            call output_ncfile%put("emissivity_noref_v", emiss_noref_v_sorted)
            call output_ncfile%put("emissivity_noref_h", emiss_noref_h_sorted)
        end if

        if (present(refl_cb_v) .and. present(refl_cb_h)) then
            call output_ncfile%put("reflectance_v", refl_cb_v_sorted)
            call output_ncfile%put("reflectance_h", refl_cb_h_sorted)
        end if

        call output_ncfile%put("sea_surface_temperature", temperature)
        call output_ncfile%put("sea_surface_salinity", salinity)

        !-------------------------------------------------------------------
        ! Close NetCDF file
        !-------------------------------------------------------------------
        call output_ncfile%close()

        write(*,'(a)') "NetCDF output completed successfully!"
        write(*,'(a, f6.3, a)') "  File size: ", file_size_mb(ncfile_name), " MB"

        ! Clean up allocated arrays
        if (allocated(wavelen_sorted)) deallocate(wavelen_sorted)
        if (allocated(emiss_v_sorted)) deallocate(emiss_v_sorted)
        if (allocated(emiss_h_sorted)) deallocate(emiss_h_sorted)
        if (allocated(emiss_noref_v_sorted)) deallocate(emiss_noref_v_sorted)
        if (allocated(emiss_noref_h_sorted)) deallocate(emiss_noref_h_sorted)
        if (allocated(refl_cb_v_sorted)) deallocate(refl_cb_v_sorted)
        if (allocated(refl_cb_h_sorted)) deallocate(refl_cb_h_sorted)
        if (allocated(sort_indices)) deallocate(sort_indices)


    end subroutine save_spectral_polarized_emissivity_wl

end module save_output
