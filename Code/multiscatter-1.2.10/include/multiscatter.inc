! multiscatter.inc - Include file for Fortran interface to multiple scattering algorithms
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

! Possible status codes returned the functions below
integer, parameter :: MS_SUCCESS = 0 ! Function executed successfully
integer, parameter :: MS_FAILURE = 1 ! A generic error occurred
integer, parameter :: MS_MEMORY_ALLOCATION_ERROR = 2 ! Error occurred allocating memory
integer, parameter :: MS_INPUT_ERROR = 3 ! Problem with the input or config data
integer, parameter :: MS_INVALID_CONTEXT = 4 ! Specified context does not exist

! The ms_set_algorithms subroutine (see below) takes three arguments,
! the context, a code indicating the small-angle algorithm and a code
! indicating the wide-angle algorithm.  Codes for small-angle
! algorithm (PVC = photon variance-covariance method):
integer, parameter :: MS_SINGLE_AND_SMALL_ANGLE_NONE= 0
integer, parameter :: MS_SINGLE_SCATTERING          = 1
integer, parameter :: MS_SMALL_ANGLE_PVC_ORIGINAL   = 2  ! Hogan (Applied Optics, 2006)
integer, parameter :: MS_SMALL_ANGLE_PVC_FAST       = 3  ! Hogan (J Atmos Sci, 2008)
integer, parameter :: MS_SMALL_ANGLE_PVC_EXPLICIT   = 4  ! Eloranta-like but PVC framework

!   Codes for wide-angle algorithm (TDTS = time-dependent two-stream method):
integer, parameter :: MS_WIDE_ANGLE_NONE   = 0
integer, parameter :: MS_WIDE_ANGLE_TDTS_FORWARD_LOBE  = 1  ! Hogan and Battaglia (J Atmos Sci, 2008)
integer, parameter :: MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE  = 2  ! Radar version of TDTS

interface
   ! FUNCTIONS TO CONFIGURE THE MULTIPLE SCATTERING CALCULATION

   ! The configuration information for a multiple scattering
   ! calculations is stored in a "context", and any number of contexts
   ! may be used by the same program.  Typically a program would
   ! obtain a new context ID with the function "ms_new_context",
   ! configure it with a set of additional subroutine calls, then call
   ! one of the multiscatter* functions to run the algorithm on a
   ! profile of input data.  The multiscatter* functions may be called
   ! again with the same context, with or without changing the
   ! configuration information (such as the altitude of the instrument
   ! which in the case of a satellite instrument may change for each
   ! profile) in between. 

   ! The first function creates a new context and return its integer
   ! ID. If wavelength is less than 100.0e-6 metres (100 microns) then
   ! the instrument is assumed to be a lidar, with a top-hat receiver
   ! type, small-angle algorithm set to MS_SMALL_ANGLE_PVC_FAST and
   ! wide-angle algorithm set to
   ! MS_WIDE_ANGLE_TDTS_FORWARD_LOBE. Otherwise, the instrument is
   ! assumed to be a radar, with a Gaussian receiver type, small-angle
   ! algorithm set to MS_SINGLE_SCATTERING and wide-angle algorithm
   ! set to MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE.
   function ms_new_context(wavelength, rho_transmitter, nfov, rho_receiver) &
        bind(C, name="ms_new_context_")
     use, intrinsic :: ISO_C_BINDING
     ! This function returns the integer ID of the context:
     integer(C_INT) :: ms_new_context
     ! Wavelength in metres:
     real(C_DOUBLE), intent(in) :: wavelength
     ! Transmitter beamwidth (1/e half-angle) in radians:
     real(C_DOUBLE), intent(in) :: rho_transmitter
     ! Number of receiver fields of view:
     integer(C_INT), intent(in) :: nfov
     ! The half-angle widths of each of receiver field of view, in
     ! radians:
     real(C_DOUBLE), intent(in), dimension(*) :: rho_receiver
   end function ms_new_context

   ! Free memory allocated to a particular context; this context can
   ! no longer be accessed (an MS_INVALID_CONTEXT code will be
   ! returned by the multiscatter* functions if you try)
   subroutine ms_free_context(context) &
        bind(C, name="ms_free_context_")
     use, intrinsic :: ISO_C_BINDING
     integer(C_INT), intent(in) :: context
   end subroutine ms_free_context

   ! Set the altitude of the instrument associated with the specified
   ! context
   subroutine ms_set_altitude(context, altitude) &
        bind(C, name="ms_set_altitude_")
     use, intrinsic :: ISO_C_BINDING
     integer(C_INT), intent(in) :: context
     ! Altitude of instrument in metres
     real(C_DOUBLE), intent(in) :: altitude
   end subroutine ms_set_altitude

   ! Set the receiver type associated with the specified context to
   ! "top hat" (usual for lidar)
   subroutine ms_set_top_hat_receiver(context) &
        bind(C, name="ms_set_top_hat_receiver_")
     use, intrinsic :: ISO_C_BINDING
     integer(C_INT), intent(in) :: context
   end subroutine ms_set_top_hat_receiver
 
   ! Set the receiver type associated with the specified context to
   ! "Gaussian" (usual for radar)
   subroutine ms_set_gaussian_receiver(context) &
        bind(C, name="ms_set_gaussian_receiver_")
     use, intrinsic :: ISO_C_BINDING
     integer(C_INT), intent(in) :: context
   end subroutine ms_set_gaussian_receiver

   ! Set the coherent backscatter enhancement; Michael Mishchenko
   ! believes that this should be 2 for spaceborne radar and 1 for
   ! lidar
   subroutine ms_set_coherent_backscatter_enhancement(context, cbe) &
        bind(C, name="ms_set_coherent_backscatter_enhancement_")
     use, intrinsic :: ISO_C_BINDING
     integer(C_INT), intent(in) :: context
     ! Coherent backscatter enhancement (1.0-2.0):
     real(C_DOUBLE), intent(in) :: cbe
   end subroutine ms_set_coherent_backscatter_enhancement

   ! Set the small-angle and wide-angle algorithms associated with the
   ! specified context
   subroutine ms_set_algorithms(context, &
        small_angle_algorithm, wide_angle_algorithm) &
        bind(C, name="ms_set_algorithms_")
     use, intrinsic :: ISO_C_BINDING
     integer(C_INT), intent(in) :: context
     ! Small-angle algorithm: should be one of
     ! MS_SINGLE_AND_SMALL_ANGLE_NONE, MS_SINGLE_SCATTERING,
     ! MS_SMALL_ANGLE_PVC_ORIGINAL, MS_SMALL_ANGLE_PVC_FAST or
     ! MS_SMALL_ANGLE_PVC_EXPLICIT (definitions above):
     integer(C_INT), intent(in) :: small_angle_algorithm
     ! Wide-angle algorithm: should be one of MS_WIDE_ANGLE_NONE,
     ! MS_WIDE_ANGLE_TDTS_FORWARD_LOBE or
     ! MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE (definitions above):
     integer(C_INT), intent(in) :: wide_angle_algorithm
   end subroutine ms_set_algorithms
 
   ! Set the small-angle algorithm to MS_SMALL_ANGLE_PVC_EXPLICIT and
   ! set the maximum number of small-angle scattering orders.  For
   ! larger numbers of scattering orders (e.g. 7), this scheme is slow
   ! but should be accurate, and is equivalent to Eloranta's method
   ! but in a PVC framework.  Note that the *_FAST and *_ORIGINAL
   ! algorithms in principle treat all orders of scattering, but in an
   ! approximate way, so may be more accurate than the *_EXPLICIT
   ! algorithm if used with an insufficient number of scattering
   ! orders.
   subroutine ms_set_explicit_max_scattering_order(context, &
        max_scattering_order) &
        bind(C, name="ms_set_explicit_max_scattering_order_")
     use, intrinsic :: ISO_C_BINDING
     integer(C_INT), intent(in) :: context
     integer(C_INT), intent(in) :: max_scattering_order
   end subroutine ms_set_explicit_max_scattering_order
    
   ! Assuming that ms_new_context was called with nfov (the number of
   ! receiver fields of view) greater than 1, this subroutine
   ! specifies that the second and later fields of view are to be
   ! treated as annular detectors, so not including the radiation from
   ! the earlier fields of view
   subroutine ms_set_annular_detectors(context) &
        bind(C, name="ms_set_annular_detectors_")
     use, intrinsic :: ISO_C_BINDING
     integer(C_INT), intent(in) :: context
   end subroutine ms_set_annular_detectors
   
   ! Only start the wide-angle calculation after a depth into the
   ! cloud when the mean-free-path is small enough compared to the
   ! field of view that this contribution to the multiple scattering
   ! return is significant.  This option speeds the calculation up in
   ! optically thin clouds, but is less accurate.
   subroutine ms_optimize_wide_angle_gates(context) &
        bind(C, name="ms_optimize_wide_angle_gates_")
     use, intrinsic :: ISO_C_BINDING
     integer(C_INT), intent(in) :: context
   end subroutine ms_optimize_wide_angle_gates

   ! Do not write messages to standard output or standard error,
   ! except for error messages
   subroutine ms_set_quiet(context) &
        bind(C, name="ms_set_quiet_")
     use, intrinsic :: ISO_C_BINDING
     integer(C_INT), intent(in) :: context
   end subroutine ms_set_quiet

   ! FUNCTIONS FOR RUNNING THE MULTIPLE SCATTERING ALGORITHM

   ! Run the algorithm using the full range of inputs, with the output
   ! backscatter containing the sum of the air and the particulate
   ! backscatter.  Note that bscat_out will need to be m*nfov in
   ! length, where nfov is the number of fields of view specified in
   ! the call to ms_new_context.
   function multiscatter(context, n, m, range, radius, &
        ext, ssa, g, ext_bscat_ratio, ext_air, ssa_air, &
        droplet_fraction, pristine_ice_fraction, &
        bscat_out) &
        bind(C, name="multiscatter_")
     use, intrinsic :: ISO_C_BINDING
     ! This function returns one of the status codes above
     integer(C_INT) :: multiscatter
     ! The "context" obtained from a call to ms_new_context
     integer(C_INT), intent(in) :: context
     ! Number of input gates
     integer(C_INT), intent(in) :: n
     ! Number of output gates (>=n)
     integer(C_INT), intent(in) :: m
     ! Height of each range gate, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: range
     ! Particle equivalent-area radius, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: radius
     ! Particle extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext
     ! Particle single scattering albedo:
     real(C_DOUBLE), intent(in), dimension(n) :: ssa
     ! Particle asymmetry factor:
     real(C_DOUBLE), intent(in), dimension(n) :: g
     ! Particle extinction-to-backscatter ratio, sr:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_bscat_ratio
     ! Air extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_air
     ! Air single scattering albedo:
     real(C_DOUBLE), intent(in), dimension(n) :: ssa_air
     ! Fraction of particle extinction from droplets (0-1):
     real(C_DOUBLE), intent(in), dimension(n) :: droplet_fraction
     ! Fraction of particle extinction  from pristine ice (0-1):
     real(C_DOUBLE), intent(in), dimension(n) :: pristine_ice_fraction
     ! Computed backscatter of the particles and air molecules
     ! combined, m-1 sr-1; this will have the dimensions of m
     ! multiplied by the number of receiver fields of view
     real(C_DOUBLE), intent(out), dimension(*) :: bscat_out
   end function multiscatter

   ! As "multiscatter", but separating out the particulate and air
   ! returns as seen by an HSRL lidar
   function multiscatter_hsrl(context, n, m, range, radius, &
        ext, ssa, g, ext_bscat_ratio, ext_air, ssa_air, &
        droplet_fraction, pristine_ice_fraction, &
        bscat_out, bscat_air_out) &
        bind(C, name="multiscatter_hsrl_")
     use, intrinsic :: ISO_C_BINDING
     ! This function returns one of the status codes above
     integer(C_INT) :: multiscatter_hsrl
     ! The "context" obtained from a call to ms_new_context
     integer(C_INT), intent(in) :: context
     ! Number of input gates
     integer(C_INT), intent(in) :: n
     ! Number of output gates (>=n)
     integer(C_INT), intent(in) :: m
     ! Height of each range gate, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: range
     ! Particle equivalent-area radius, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: radius
     ! Particle extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext
     ! Particle single scattering albedo:
     real(C_DOUBLE), intent(in), dimension(n) :: ssa
     ! Particle asymmetry factor:
     real(C_DOUBLE), intent(in), dimension(n) :: g
     ! Particle extinction-to-backscatter ratio, sr:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_bscat_ratio
     ! Air extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_air
     ! Air single scattering albedo:
     real(C_DOUBLE), intent(in), dimension(n) :: ssa_air
     ! Fraction of particle extinction from droplets (0-1):
     real(C_DOUBLE), intent(in), dimension(n) :: droplet_fraction
     ! Fraction of particle extinction  from pristine ice (0-1):
     real(C_DOUBLE), intent(in), dimension(n) :: pristine_ice_fraction
     ! Computed backscatter of the particles only, m-1 sr-1; this will
     ! have the dimensions of m multiplied by the number of receiver
     ! fields of view
     real(C_DOUBLE), intent(out), dimension(*) :: bscat_out
     ! Computed backscatter of the air molecules only, m-1 sr-1; this
     ! will have the dimensions of m multiplied by the number of
     ! receiver fields of view
     real(C_DOUBLE), intent(out), dimension(*) :: bscat_air_out
   end function multiscatter_hsrl
 
   ! As "multiscatter", but without the inputs necessary for the
   ! wide-angle scattering calculation; thus, even if a wide-angle
   ! algorithm has been specified, this part of the algorithm will not
   ! be called. This function would typically be used for lidar
   ! observations of ice clouds where the wide-angle scattering is
   ! typically very small, and saves having to call "multiscatter"
   ! with dummy arguments for the unused variables.
   function multiscatter_small_angle_only(context, n, range, radius, &
        ext, ext_bscat_ratio, ext_air, &
        droplet_fraction, pristine_ice_fraction, &
        bscat_out) &
        bind(C, name="multiscatter_small_angle_only_")
     use, intrinsic :: ISO_C_BINDING
     ! This function returns one of the status codes above
     integer(C_INT) :: multiscatter_small_angle_only
     ! The "context" obtained from a call to ms_new_context
     integer(C_INT), intent(in) :: context
     ! Number of input and output gates
     integer(C_INT), intent(in) :: n
     ! Height of each range gate, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: range
     ! Particle equivalent-area radius, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: radius
     ! Particle extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext
     ! Particle extinction-to-backscatter ratio, sr:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_bscat_ratio
     ! Air extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_air
     ! Fraction of particle extinction from droplets (0-1):
     real(C_DOUBLE), intent(in), dimension(n) :: droplet_fraction
     ! Fraction of particle extinction  from pristine ice (0-1):
     real(C_DOUBLE), intent(in), dimension(n) :: pristine_ice_fraction
     ! Computed backscatter of the particles and air molecules
     ! combined, m-1 sr-1; this will have the dimensions of m
     ! multiplied by the number of receiver fields of view
     real(C_DOUBLE), intent(out), dimension(*) :: bscat_out
   end function multiscatter_small_angle_only

   ! As multiscatter_small_angle but with separate outputs for the
   ! particulate and the air backscatter
   function multiscatter_small_angle_only_hsrl(context, n, range, radius, &
        ext, ext_bscat_ratio, ext_air, &
        droplet_fraction, pristine_ice_fraction, &
        bscat_out, bscat_air_out) &
        bind(C, name="multiscatter_small_angle_only_hsrl_")
     use, intrinsic :: ISO_C_BINDING
     ! This function returns one of the status codes above
     integer(C_INT) :: multiscatter_small_angle_only_hsrl
     ! The "context" obtained from a call to ms_new_context
     integer(C_INT), intent(in) :: context
     ! Number of input and output gates
     integer(C_INT), intent(in) :: n
     ! Height of each range gate, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: range
     ! Particle equivalent-area radius, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: radius
     ! Particle extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext
     ! Particle extinction-to-backscatter ratio, sr:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_bscat_ratio
     ! Air extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_air
     ! Fraction of particle extinction from droplets (0-1):
     real(C_DOUBLE), intent(in), dimension(n) :: droplet_fraction
     ! Fraction of particle extinction  from pristine ice (0-1):
     real(C_DOUBLE), intent(in), dimension(n) :: pristine_ice_fraction
     ! Computed backscatter of the particles and air molecules
     ! combined, m-1 sr-1; this will have the dimensions of m
     ! multiplied by the number of receiver fields of view
     real(C_DOUBLE), intent(out), dimension(*) :: bscat_out
     ! Computed backscatter of the air molecules only, m-1 sr-1; this
     ! will have the dimensions of m multiplied by the number of
     ! receiver fields of view
     real(C_DOUBLE), intent(out), dimension(*) :: bscat_air_out
   end function multiscatter_small_angle_only_hsrl
 
   ! A simple interface that includes only small-angle scattering and
   ! assumes that the near-backscatter phase function is isotropic, so
   ! the parameterizations embodied in droplet_fraction and
   ! pristine_ice_fraction are not required
   function multiscatter_simple(context, n, range, radius, &
        ext, ext_bscat_ratio, ext_air, &
        bscat_out) &
        bind(C, name="multiscatter_simple_")
     use, intrinsic :: ISO_C_BINDING
     ! This function returns one of the status codes above
     integer(C_INT) :: multiscatter_simple
     ! The "context" obtained from a call to ms_new_context
     integer(C_INT), intent(in) :: context
     ! Number of input and output gates
     integer(C_INT), intent(in) :: n
     ! Height of each range gate, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: range
     ! Particle equivalent-area radius, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: radius
     ! Particle extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext
     ! Particle extinction-to-backscatter ratio, sr:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_bscat_ratio
     ! Air extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_air
     ! Computed backscatter of the particles only, m-1 sr-1; this will
     ! have the dimensions of m multiplied by the number of receiver
     ! fields of view
     real(C_DOUBLE), intent(out), dimension(*) :: bscat_out
   end function multiscatter_simple
   
   ! The equivalent adjoint code for the "multiscatter" function, used
   ! in quasi-Newton retrieval schemes
   function multiscatter_AD(context, n, m, range, radius, &
        ext, ssa, g, ext_bscat_ratio, ext_air, ssa_air, &
        droplet_fraction, pristine_ice_fraction, &
        bscat_out, bscat_AD, radius_AD, &
        ext_AD, ssa_AD, g_AD, ext_bscat_ratio_AD, ext_air_AD, &
        droplet_fraction_AD, pristine_ice_fraction_AD) &
        bind(C, name="multiscatter_AD_")
     use, intrinsic :: ISO_C_BINDING
     ! This function returns one of the status codes above
     integer(C_INT) :: multiscatter_AD
     ! The "context" obtained from a call to ms_new_context
     integer(C_INT), intent(in) :: context
     ! Number of input gates
     integer(C_INT), intent(in) :: n
     ! Number of output gates (>=n)
     integer(C_INT), intent(in) :: m
     ! Height of each range gate, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: range
     ! Particle equivalent-area radius, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: radius
     ! Particle extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext
     ! Particle single scattering albedo:
     real(C_DOUBLE), intent(in), dimension(n) :: ssa
     ! Particle asymmetry factor:
     real(C_DOUBLE), intent(in), dimension(n) :: g
     ! Particle extinction-to-backscatter ratio, sr:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_bscat_ratio
     ! Air extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_air
     ! Air single scattering albedo:
     real(C_DOUBLE), intent(in), dimension(n) :: ssa_air
     ! Fraction of particle extinction from droplets (0-1):
     real(C_DOUBLE), intent(in), dimension(n) :: droplet_fraction
     ! Fraction of particle extinction  from pristine ice (0-1):
     real(C_DOUBLE), intent(in), dimension(n) :: pristine_ice_fraction
     ! Computed backscatter of the particles and air molecules
     ! combined, m-1 sr-1; this will have the dimensions of m
     ! multiplied by the number of receiver fields of view
     real(C_DOUBLE), intent(out), dimension(*) :: bscat_out
     ! Adjoint input (same dimension as bscat_out)
     real(C_DOUBLE), intent(in), dimension(*) :: bscat_AD
     ! Adjoint outputs (modifies the input variable)
     real(C_DOUBLE), intent(inout), dimension(n) :: radius_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: ext_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: ssa_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: g_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: ext_bscat_ratio_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: ext_air_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: droplet_fraction_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: pristine_ice_fraction_AD
   end function multiscatter_AD

   ! The equivalent adjoint code for the "multiscatter_hsrl" function, used
   ! in quasi-Newton retrieval schemes
   function multiscatter_hsrl_AD(context, n, m, range, radius, &
        ext, ssa, g, ext_bscat_ratio, ext_air, ssa_air, &
        droplet_fraction, pristine_ice_fraction, &
        bscat_out, bscat_air_out, bscat_AD, bscat_air_AD, radius_AD, &
        ext_AD, ssa_AD, g_AD, ext_bscat_ratio_AD, ext_air_AD, &
        droplet_fraction_AD, pristine_ice_fraction_AD) &
        bind(C, name="multiscatter_hsrl_AD_")
     use, intrinsic :: ISO_C_BINDING
     ! This function returns one of the status codes above
     integer(C_INT) :: multiscatter_hsrl_AD
     ! The "context" obtained from a call to ms_new_context
     integer(C_INT), intent(in) :: context
     ! Number of input gates
     integer(C_INT), intent(in) :: n
     ! Number of output gates (>=n)
     integer(C_INT), intent(in) :: m
     ! Height of each range gate, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: range
     ! Particle equivalent-area radius, metres:
     real(C_DOUBLE), intent(in), dimension(n) :: radius
     ! Particle extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext
     ! Particle single scattering albedo:
     real(C_DOUBLE), intent(in), dimension(n) :: ssa
     ! Particle asymmetry factor:
     real(C_DOUBLE), intent(in), dimension(n) :: g
     ! Particle extinction-to-backscatter ratio, sr:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_bscat_ratio
     ! Air extinction coefficient, m-1:
     real(C_DOUBLE), intent(in), dimension(n) :: ext_air
     ! Air single scattering albedo:
     real(C_DOUBLE), intent(in), dimension(n) :: ssa_air
     ! Fraction of particle extinction from droplets (0-1):
     real(C_DOUBLE), intent(in), dimension(n) :: droplet_fraction
     ! Fraction of particle extinction  from pristine ice (0-1):
     real(C_DOUBLE), intent(in), dimension(n) :: pristine_ice_fraction
     ! Computed backscatter of the particles, m-1 sr-1; this will have
     ! the dimensions of m multiplied by the number of receiver fields
     ! of view
     real(C_DOUBLE), intent(out), dimension(*) :: bscat_out
     ! Computed backscatter of the air molecules, m-1 sr-1; this will
     ! have the dimensions of m multiplied by the number of receiver
     ! fields of view 
     real(C_DOUBLE), intent(out), dimension(*) :: bscat_air_out
     ! Adjoint inputs (same dimensions as bscat_out)
     real(C_DOUBLE), intent(in), dimension(*) :: bscat_AD
     real(C_DOUBLE), intent(in), dimension(*) :: bscat_air_AD
     ! Adjoint outputs (modifies the input array)
     real(C_DOUBLE), intent(inout), dimension(n) :: radius_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: ext_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: ssa_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: g_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: ext_bscat_ratio_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: ext_air_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: droplet_fraction_AD
     real(C_DOUBLE), intent(inout), dimension(n) :: pristine_ice_fraction_AD
   end function multiscatter_hsrl_AD

end interface
