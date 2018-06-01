! multiscatterf_ascii.f90 -- Single-profile fortran interface for
! lidar multiple scattering algorithm
!
! Copyright (C) 2010-2011 Robin Hogan <r.j.hogan@reading.ac.uk>
!
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 2 of the License, or (at your option) any later version.
!
!  This library is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!   
!  You should have received a copy of the GNU Lesser General Public
!  License along with this library; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!
! Demo of Fortran interface to lidar multiple scattering algorithm;
! this is like the C version (multiscatter_ascii.c) except that there
! are no command-line options to change the behaviour, and it is
! required that all input variables are provided
!  
! Input data (standard input)
!   First line:
!     num_gates wavelength(m) lidar_alt(m) half_angle_beam_divergence(rad) half_angle_receiver_fov(rad)
!   Subsequent lines:
!     range(m) ext(m-1) radius(m) ext_bscat_ratio(sr) ext_air(m) ssa_air(m) droplet_fraction(m) pristine_ice_fraction(m)
!
! Output data (standard output)
!     igate range(m) ext(m-1) radius(m) bscat_total(m-1sr-1)

program multiscatterf_ascii
  implicit none
  
  ! Fortran 2003 interface is defined here:
  include 'multiscatter.inc'

  ! If your compiler doesn't support the Fortran 2003 interface, include this instead
!  include 'multiscatter77.inc'
  
  integer :: n              ! Number of range-gates in profile
  integer :: i              ! Counter
  ! Multiscatter context ID
  integer :: id;
  ! The status (error code) of an algorithm calculation
  integer :: status
  logical :: is_radar = .false. ! Use wavelength to guess if is a radar or lidar

  integer, parameter :: MAX_LEN = 1000
  real*8 :: wavelength      ! Wavelength in m
  real*8 :: rho_transmitter ! Transmitter 1/e half width in radians
  real*8 :: rho_receiver(1) ! Receiver half width in radians
  real*8 :: altitude        ! Instrument altitude in m
  real*8 :: range(MAX_LEN)  ! Range to each range-gate, m
  real*8 :: radius(MAX_LEN) ! Equivalent-area radius, m
  real*8 :: ext(MAX_LEN)    ! Particle extinction coefficient, m-1
  real*8 :: ssa(MAX_LEN)    ! Particle single scatteing albedo
  real*8 :: g(MAX_LEN)      ! Particle asymmetry factor
  real*8 :: ext_bscat_ratio(MAX_LEN) ! Particle extinction-backscatter ratio, sr
  real*8 :: ext_air(MAX_LEN)! Air extinction coefficient, m-1
  real*8 :: ssa_air(MAX_LEN)! Air single scattering albedo
  real*8 :: droplet_fraction(MAX_LEN)! Fraction of particulate extinction due to droplets
  real*8 :: pristine_ice_fraction(MAX_LEN)! Frac of part. ext. due to pristine ice
  real*8 :: bscat(MAX_LEN)  ! Backscatter output from algorithm, m-1 sr-1

  character (len=300) :: line
 
  ! Skip over lines beginning with "#"
  read_comments: do
     read(*, '(A)') line
     if (line(1:1) /= '#') exit read_comments
  end do read_comments
  
  ! First line containing instrument properties must be read from a string
  read(line,*) n, wavelength, altitude, rho_transmitter, rho_receiver(1)

  if (n > MAX_LEN) then
     write(*,*) '*** Error: number of range gates (', n, ') exceeds limit (', MAX_LEN, ') ***'
     stop
  end if

  ! If this is a radar then we will choose different algorithms
  if (wavelength > 100.0e-6) then
     is_radar = .true.
  end if

  
  ! Load subsequent lines
  do i = 1,n
     read(*,*) range(i), ext(i), radius(i), &
          ext_bscat_ratio(i), ext_air(i), ssa(i), g(i), ssa_air(i), &
          droplet_fraction(i), pristine_ice_fraction(i)
  end do
  

  !  Get a "context" for the calculation, into which we will put
  !  information about the problem to be solved
  id = ms_new_context(wavelength, rho_transmitter, &
       1, rho_receiver)

  !  Set the altitude and the algorithm to use
  call ms_set_altitude(id, altitude)
  if (is_radar) then
     call ms_set_algorithms(id, &
          MS_SINGLE_SCATTERING, MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE)
  else 
     call ms_set_algorithms(id, &
          MS_SMALL_ANGLE_PVC_FAST, MS_WIDE_ANGLE_TDTS_FORWARD_LOBE)
  end if


  ! We may wish to only perform the wide-angle calculation if the
  ! mean-free-path is small enough, in which case uncomment this line:
  !   call ms_optimize_wide_angle_gates(id);

  !  Call the algorithm
  status = multiscatter(id, &
            n, n, range, radius, ext, ssa, g, ext_bscat_ratio, &
            ext_air, ssa_air, droplet_fraction, pristine_ice_fraction, &
            bscat)
  if (status /= MS_SUCCESS) then
     write(*,*) '*** An error occurred ***'
     stop
  end if
  

  !  Write the results to standard output
  do i = 1,n
     write(*,'(I4, 7E11.4)') i, range(i), ext(i), radius(i), &
          bscat(i)
  end do
  
  !  Free allocated data associated with context "id"
  call ms_free_context(id)

  stop
end program multiscatterf_ascii

