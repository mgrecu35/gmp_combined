/* multiscatter_small_angle_ADonly.c -- Approximate adjoint of a
   small-angle multiple scattering calculation

   Copyright (C) 2010 Robin Hogan <r.j.hogan@reading.ac.uk> 

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

#include <math.h>
#include "multiscatter.h"

int
multiscatter_small_angle_ADonly(int n, ms_instrument instrument,
				ms_surface surface, ms_real *range,
				ms_real *ext, ms_real *ext_bscat_ratio,
				ms_real *ext_air, ms_real *ssa_air,
				ms_real *bscat, ms_real *bscat_air,
				ms_real *trans_enhancement,
				ms_real *anisotropic_factor,
				/* Adjoint terms (in) */
				ms_real *bscat_AD, ms_real *bscat_air_AD,
				/* Adjoint terms (out) */
				ms_real *ext_AD, ms_real *ext_bscat_ratio_AD)
{
  ms_real platt_factor[n]; /* Platt's multiple-scattering factor */
  ms_real eff_transmittance[n]; /* Effective transmittance to end of
				   layer including multiple scattering
				   enhancement */
  ms_real optical_depth = 0.0; /* Two-way optical depth to start of current layer */
  ms_real two_drange = 2.0*fabs(range[1]-range[0]);
  ms_real bscat_ext_ratio_air = 3.0 / (8.0 * MS_PI);
  int i;
  
  /* Estimate the equivalent Platt multiple-scattering factor versus
     range */
  for (i = 0; i < n; i++) {
    ms_real od; /* Two-way optical depth of layer i */
    if (ext_air) {
      od = two_drange*(ext[i]+ext_air[i]);
    }
    else {
      od = two_drange*ext[i];
    }

    optical_depth += od;
    
    if (optical_depth > 0.0) {
      platt_factor[i] = 1.0 - log(1.0+trans_enhancement[i]) / optical_depth;
    }
    else {
      platt_factor[i] = 1.0;
    }

    eff_transmittance[i] = exp(-optical_depth) * (1.0+trans_enhancement[i]);
  }

  /* *** ADJOINT CALCULATION *** */
  if (ext_air && ssa_air) {
    /* Scattering atmospheric gases are present */
    ms_real sum = 0.0;
    for (i = n-1; i >= 0; i--) {
      ms_real B = 0.0;
      ms_real Bair = 0.0;
      ms_real eff_trans_prev, factor;
      if (i > 0) {
	eff_trans_prev = eff_transmittance[i-1];
      }
      else {
	eff_trans_prev = 1.0;
      }
      factor = two_drange*platt_factor[i]*eff_transmittance[i]
	/ /* (1.0+trans_enhancement[i]) *         */
	(eff_trans_prev-eff_transmittance[i])
	- 1.0/(ext[i] + ext_air[i]);


      if (bscat_air) {
	/* Particulate and air backscatters are being separated */
	if (ext[i] > 0.0) {
	  B = bscat[i] * (1.0/ext[i] + factor);
	}
	else {
	  B = eff_trans_prev/ext_bscat_ratio[i];
	}

	if (ext_air[i] > 0.0) {
	  /*	  Bair = bscat_air[i] * factor;   */
	  Bair = (ext_air[i]*ssa_air[i]*bscat_ext_ratio_air
		  *eff_transmittance[i] - bscat_air[i])
	    /(ext[i]+ext_air[i]);
	}

	ext_AD[i] += B*bscat_AD[i] + Bair*bscat_air_AD[i] - sum;

	sum += two_drange*platt_factor[i]*(bscat[i]*bscat_AD[i]
					 + bscat_air[i]*bscat_air_AD[i]);
	ext_bscat_ratio_AD[i] -= bscat_AD[i]*bscat[i]/ext_bscat_ratio[i];
      }
      else {
	/* "bscat" and "bscat_AD" refer to the total backscatter */
	if (ext[i] > 0.0) {
	  B = bscat[i] * (1.0/(ext[i] + bscat_ext_ratio_air*ssa_air[i]
			      *ext_air[i]*ext_bscat_ratio[i])
			  + factor);
	}

	ext_AD[i] += B*bscat_AD[i] - sum;
	sum += two_drange*platt_factor[i]*(bscat[i]*bscat_AD[i]);
	if (ext[i] > 0.0) {
	  /* This neglects anisotropic phase function effects, but
	     simply weights the contribution by their unattenuated
	     backscatters: */
	  ext_bscat_ratio_AD[i] -= bscat_AD[i]*bscat[i]
	    /(ext_bscat_ratio[i]*(1.0 + ext_bscat_ratio[i]
			  *bscat_ext_ratio_air*ssa_air[i]*ext_air[i]/ext[i]));
	}	
      }
    }
  }
  else { /* THIS BIT OF CODE HAS NOT BEEN CHECKED */
    /* No atmospheric gases to consider for scattering, although they
       might attenuate */ 
    ms_real sum = 0.0;
    for (i = n-1; i >= 0; i--) {
      if (ext[i] > 0.0) {
	ms_real B;

	if (ext_air) {
	  ms_real eff_trans_prev, factor;
	  if (i > 0) {
	    eff_trans_prev = eff_transmittance[i-1];
	  }
	  else {
	    eff_trans_prev = 1.0;
	  }

      /*
	  ms_real factor = two_drange*eff_transmittance[i]
	    / ((1.0+trans_enhancement[i])
	       *(eff_trans_prev-eff_transmittance[i]))
	    - 1.0/(ext[i] + ext_air[i]);
      */
	  factor = two_drange*platt_factor[i]*eff_transmittance[i]
	    / (eff_trans_prev-eff_transmittance[i])
	    - 1.0/(ext[i] + ext_air[i]);
	  B = bscat[i] * (1.0/ext[i] + factor);
	}
	else {
	  B = bscat[i] / ext[i];
	}
	ext_AD[i] += B*bscat_AD[i] - sum;
	sum += two_drange*platt_factor[i]*bscat[i]*bscat_AD[i];
	ext_bscat_ratio_AD[i] -= bscat_AD[i]*bscat[i]/ext_bscat_ratio[i];
      }
    }
  }

  return MS_SUCCESS;
}
