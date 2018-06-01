/* multiscatter.c -- Adjoint of interface to multiple scattering algorithms

   Copyright (C) 2006-2010 Robin Hogan <r.j.hogan@reading.ac.uk> 

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


   The algorithm is implemented in ANSI C, but a Fortran interface is
   provided.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* For details of how to call the multiple scattering algorithms, see
   the comments in this header file: */
#include "ms.h"


/* Perform the multiple scattering calculation, using which ever
   combination of algorithms is appropriate depending on multiple
   scattering regime */
int
multiscatter_AD(
    /* Input data */
    int n,                    /* Number of input gates */
    int m,                    /* Number of output gates (>= n) */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,    /* Cloud/aerosol equivalent radius, microns */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./backscatter ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft., m-1 (NULL for vacuum) */
    const ms_real *ssa_air,   /* Air single-scatter albedo */
    const ms_real *droplet_fraction,     /* Fraction of extinction from droplets */
    const ms_real *pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */
    ms_real *bscat_air_out,   /* Measured backscatter of air, m-1 sr-1 */
    /* Adjoint inputs */
    const ms_real *bscat_AD,
    const ms_real *bscat_air_AD,
    /* Adjoint outputs */
    ms_real *radius_AD,
    ms_real *ext_AD,
    ms_real *ssa_AD,
    ms_real *g_AD,
    ms_real *ext_bscat_ratio_AD,
    ms_real *ext_air_AD,
    ms_real *droplet_fraction_AD,
    ms_real *pristine_ice_fraction_AD)
{
  ms_instrument instrument_1fov = instrument;
  ms_real ss_multiplier = 1.0;
  int status = MS_SUCCESS;
  int i;
  int ifov;
  int require_bscat_air = (bscat_air_out != NULL);
  int nan_input = 0;

  if (config->options & MS_CHECK_FOR_NAN
      && ms_isnan(n, ext_AD)) {
    fprintf(stderr, "Warning: ext_AD input to multiscatter_AD() contains NaNs\n");
    nan_input = 1;
  }
  //  fprintf(stderr, "!!!");
    for (i = 0; i < n; i++) {
      //      fprintf(stderr, " %g", ext_AD[i]);
      ext_AD[i] = 0.0;
    }
    //  fprintf(stderr, "\n");

  /* This variable contains the same instrument settings, except for
     having only one field-of-view; it is modified when we loop over
     the fields-of-view */
  instrument_1fov.nfov = 1;

  /* Work out single-scattering factor: this ensures that the first
     field-of-view is properly calibrated and that the backscatter in
     the subsequent fields-of-view are calibrated to the first such
     that the backscatter is proportional to the energy received. */
  if (instrument.receiver_type == MS_TOP_HAT) {
    ss_multiplier
      = 1.0/(1.0-exp(-instrument.rho_receiver[0]*instrument.rho_receiver[0]/
		     (instrument.rho_transmitter*instrument.rho_transmitter)));
  }
  else {
    ss_multiplier
      = 1.0 + instrument.rho_transmitter*instrument.rho_transmitter
      /(instrument.rho_receiver[0]*instrument.rho_receiver[0]);
  }
  config->ss_multiplier = ss_multiplier;

  if (!(config->options & MS_QUIET) && ss_multiplier != 1.0) {
    fprintf(stderr, "Factor to ensure calibration of first FOV: %g\n",
	    ss_multiplier);
  }


  /* FIRST DO SMALL-ANGLE SCATTERING CALCULATION */

  if (config->small_angle_algorithm != MS_SINGLE_AND_SMALL_ANGLE_NONE) {
    /* In case of annular detectors we may need to manipulate the
       adjoints before passing them into the functions - here we
       allocate the storage for the manipulated vectors */
    ms_real bscat_annular_AD_data[m];
    ms_real bscat_air_annular_AD_data[m];

    /* For small-angle scattering we loop over the fields-of-view */
    for (ifov = 0; ifov < instrument.nfov; ifov++) {
      /* Set the receiver field of view appropriately */
      const ms_real* bscat_annular_AD = bscat_AD + ifov*m;
      const ms_real* bscat_air_annular_AD 
	= bscat_air_AD + require_bscat_air*ifov*m;
      ms_real rho_receiver = instrument.rho_receiver[ifov];
      ms_real fov_ss_multiplier = 1;

      instrument_1fov.rho_receiver = &rho_receiver;

      if (ifov > 0) {
	if (instrument.receiver_type == MS_TOP_HAT) {
	  fov_ss_multiplier
	    = ss_multiplier*(1.0-exp(-instrument.rho_receiver[ifov]
				     *instrument.rho_receiver[ifov]/
				     (instrument.rho_transmitter
				      *instrument.rho_transmitter)));
	}
	else {
	  fov_ss_multiplier = ss_multiplier
	    / (1.0 + instrument.rho_transmitter*instrument.rho_transmitter
	       /(instrument.rho_receiver[ifov]*instrument.rho_receiver[ifov]));
	}
      }


      /* We need to manipulate adjoints if we are using annular
	 detectors AND we have more than one field of view AND it is
	 not the final field of view */
      if (config->options & MS_ANNULAR_DETECTORS) {
	if (ifov < instrument.nfov-1) {
	  int i;
	  for (i = 0; i < m; i++) {
	    bscat_annular_AD_data[i] = 1.0/*fov_ss_multiplier*/
	      * (bscat_annular_AD[i] - bscat_annular_AD[i+m]);
	  }
	  bscat_annular_AD = bscat_annular_AD_data;
	  if (require_bscat_air) {
	    for (i = 0; i < m; i++) {
	      bscat_air_annular_AD_data[i] = 1.0/*fov_ss_multiplier*/
		* (bscat_air_annular_AD[i] - bscat_air_annular_AD[i+m]);
	    }
	    bscat_air_annular_AD = bscat_air_annular_AD_data;
	  }
	}
	/*
	else {
	  int i;
	  for (i = 0; i < m; i++) {
	    bscat_annular_AD_data[i] = fov_ss_multiplier
	      * bscat_annular_AD[i];
	  }
	  bscat_annular_AD = bscat_annular_AD_data;
	  if (require_bscat_air) {
	    for (i = 0; i < m; i++) {
	      bscat_air_annular_AD_data[i] = fov_ss_multiplier
		* bscat_air_annular_AD[i];
	    }
	    bscat_air_annular_AD = bscat_air_annular_AD_data;
	  }
	}
	*/

      }

      else if (fov_ss_multiplier != 1.0
	       && ifov > 0) {
	int i;
	for (i = 0; i < m; i++) {
	  bscat_annular_AD_data[i] = fov_ss_multiplier
	    * bscat_annular_AD[i];
	}
	bscat_annular_AD = bscat_annular_AD_data;
	if (require_bscat_air) {
	  for (i = 0; i < m; i++) {
	    bscat_air_annular_AD_data[i] = fov_ss_multiplier
	      * bscat_air_annular_AD[i];
	  }
	  bscat_air_annular_AD = bscat_air_annular_AD_data;
	}
      }

      if (config->small_angle_algorithm == MS_SINGLE_SCATTERING) {
	/* If no forward lobe (radar case) then we have no small-angle
	   multiple scattering contribution and we simply calculate
	   the single scattering here, with wide-angle scattering
	   later. */
	status
	  = ms_singlescatter_AD(n, instrument_1fov, surface,
				range, ext, ext_bscat_ratio,
				ext_air,
				bscat_out + ifov*m,
				bscat_air_out + ifov*m*require_bscat_air,
				bscat_annular_AD,
				bscat_air_annular_AD,
				ext_AD, ext_bscat_ratio_AD);
      }
      else {
	if (config->small_angle_algorithm != MS_SMALL_ANGLE_PVC_EXPLICIT) {
	  /* Perform small-angle multiple-scattering calculation using an
	     efficient algorithm */
	  if (config->small_angle_algorithm == MS_SMALL_ANGLE_PVC_FAST_LAG) {
	    fprintf(stderr, "Error: adjoint not available with small-angle lag calculation\n");
	    return MS_FAILURE;
	  }
	  else if (config->small_angle_algorithm == MS_SMALL_ANGLE_PVC_FAST) {
	    status = ms_fast_small_angle_AD(n, 
			 instrument_1fov, surface, 
			 range, radius, ext, ext_bscat_ratio, ext_air,
			 droplet_fraction, pristine_ice_fraction,
			 bscat_out + ifov*m,
			 bscat_air_out + ifov*m*require_bscat_air,
		         bscat_AD, bscat_air_AD, radius_AD, 
			 ext_AD, ext_bscat_ratio_AD, ext_air_AD,
			 droplet_fraction_AD, pristine_ice_fraction_AD);
	  }
	  else {
	    fprintf(stderr, "Error: adjoint not available for O(N^2) small-angle calculation (Hogan 2006): use the \"-fast-sa\" option\n");
	    return MS_FAILURE;
	  }
	}
	else {
	  fprintf(stderr, "Error: adjoint not available for explicit small-angle calculation: use the \"-fast-sa\" option\n");
	  return MS_FAILURE;
	}
	if (status != MS_SUCCESS) {
	  return status;
	}

      } /* End of if single scattering else small-angle scattering */

      /* To ensure the correct calibration, we may need to scale the
	 wider fields of view - here we multiply by the multiplier for
	 the first field-of-view, and divide by the multiplier
	 calculated for the current field-of-view.  At the wide-angle
	 scattering stage, we simply multiply by the multiplier for
	 the first field-of-view.  This ensures that for a given
	 field-of-view, the narrow-angle and wide-angle contributions
	 are correctly scaled relative to each other (normally one
	 would only multiply the wide-angle component by the
	 multiplier for the current FOV), and relative to the first
	 field-of-view */
      if (ifov > 0) {
	ms_real* bscat_cur = bscat_out+ifov*m;
	ms_real* bscat_air_cur = bscat_air_out+require_bscat_air*ifov*m;
	for (i = 0; i < n; i++, bscat_cur++) {
	  (*bscat_cur) *= fov_ss_multiplier;
	}
	if (bscat_air_cur) {
	  for (i = 0; i < n; i++, bscat_air_cur++) {
	    (*bscat_air_cur) *= fov_ss_multiplier;
	  }
	}
	
	if (!(config->options & MS_QUIET) && ss_multiplier != 1.0) {
	  fprintf(stderr, "   Equivalent factor for FOV %d: %g\n",
		  ifov+1, ss_multiplier/fov_ss_multiplier);
	}
      }
    } /* Loop over fields-of-view */

    if (config->options & MS_ANNULAR_DETECTORS) {
      /* Need to subtract backscatter contributions from each other */
      for (ifov = instrument.nfov-1; ifov > 0; ifov--) {
	int i;
	int i_inner = (ifov-1)*m;
	int i_outer = ifov*m;
	for (i = 0; i < m; i++, i_inner++, i_outer++) {
	  bscat_out[i_outer] -= bscat_out[i_inner];
	}
	if (bscat_air_out) {
	  i_inner = (ifov-1)*m;
	  i_outer = ifov*m;
	  for (i = 0; i < m; i++, i_inner++, i_outer++) {
	    bscat_air_out[i_outer] -= bscat_air_out[i_inner];
	  }
	}
      }
    }
  } /* If perform small-angle calculation */

  /* Return if status is non-zero (indicating that an error
     occurred) or if ssa or g are not set (indicating that the
     wide-angle calculation is not to be performed). */
  if (status != MS_SUCCESS
      || config->wide_angle_algorithm == MS_WIDE_ANGLE_NONE) {
    return status;
  }
  if ( !ssa || !g ) {
    fprintf(stderr, "Error: tried to call wide-angle part of calculation without single scattering albedo or asymmetry factor being set\n");
    return MS_FAILURE;
  }

  if (config->small_angle_algorithm == MS_SINGLE_AND_SMALL_ANGLE_NONE) {
    /* Reset backscatter to zero, since wide-angle calculation
       increments existing values */
    int i;
    for (i = 0; i < m*instrument.nfov; i++) {
      bscat_out[i] = 0.0;
    }
  }

  /* NOW DO WIDE-ANGLE SCATTERING CALCULATION WITH ADJOINT */
  status = ms_wide_angle_AD(n, m,
	    config, instrument, surface, range, radius,
	    ext, ssa, g, ext_air, ssa_air,
	    bscat_out,
	    bscat_AD, 
            radius_AD, ext_AD, ssa_AD, g_AD);

  if (status != MS_SUCCESS) {
    return status;
  }

  if (config->options & MS_CHECK_FOR_NAN
      && !nan_input
      && ms_isnan(n, ext_AD)) {
    fprintf(stderr, "Warning: algorithm produced NaNs in ext_AD\n");
    ms_dump_input(NULL, n, config, instrument, surface, range, radius,
		  ext, ssa, g, ext_bscat_ratio, ext_air, ssa_air,
		  droplet_fraction, pristine_ice_fraction,
		  bscat_AD, bscat_air_AD);
  }
  return MS_SUCCESS;
}

