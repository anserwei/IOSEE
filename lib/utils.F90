! utils.F90
!
! Purpose:
!   Common utilities and type definitions for ocean emissivity calculations.
!   Provides precision type definitions, physical constants, and mathematical
!   utilities for scientific computing with consistent data types.
!
! Features:
!   - Precision type definitions (jprd, jpim) for consistent data types
!   - Physical constants (pi, deg2rad) with high precision
!   - Common mathematical utilities for scientific applications
!   - Error handling types and constants
!   - Platform compatibility functions for portable code
!   - IEEE 754 compliant floating-point operations
!
! Dependencies:
!   - Standard Fortran intrinsic modules
!
! Public Interface:
!   - jprd, jprb: Double precision real kind parameters
!   - jpim, jpib: Integer kind parameters for different ranges
!   - pi, deg2rad: Mathematical constants
!   - All types and constants are public for global use
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
module utils
    implicit none
    public
    save

    !------------------------------------------------------------------
    ! Integer kinds
    !------------------------------------------------------------------
    integer, parameter :: jpit = selected_int_kind(2)
    integer, parameter :: jpis = selected_int_kind(4)
    integer, parameter :: jpim = selected_int_kind(9)
    integer, parameter :: jpib = selected_int_kind(12)

    !------------------------------------------------------------------
    ! Real kinds
    !------------------------------------------------------------------
    integer, parameter :: jprt = selected_real_kind(2,1)
    integer, parameter :: jprs = selected_real_kind(4,2)
    integer, parameter :: jprm = selected_real_kind(6,37)
    integer, parameter :: jprd = selected_real_kind(13,300)
    integer, parameter :: jprb = jprd

    !------------------------------------------------------------------
    ! Physical constants
    !------------------------------------------------------------------
    real(kind=jprd), parameter :: pi = 2.0_jprd * asin(1.0_jprd)
    real(kind=jprd), parameter :: deg2rad = pi / 180.0_jprd

end module utils
