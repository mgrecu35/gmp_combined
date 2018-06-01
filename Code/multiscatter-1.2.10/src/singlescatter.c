/* singlescatter.c -- Single scattering calculation

   Copyright (C) 2009 Robin Hogan <r.j.hogan@reading.ac.uk> 

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

*/

#include <math.h>

#include "ms.h"

/* Get the range gate spacing for gate i given n gates with ranges
   centred at "range". The mid-points are assumed to lie half-way
   between the centres, and the end points are treated as one would
   expect. */
ms_real
ms_get_drange(int i, int n, const ms_real *range)
{
  if (i <= 0) {
    return fabs(range[1]-range[0]);
  }
  else if (i >= n-1) {
    return fabs(range[n-1]-range[n-2]);
  }
  else {
    return fabs(range[i+1]-range[i-1])*0.5;
  }
}

/* Get the mid-point between gates i-1 and i. */
ms_real
ms_get_midpoint(int i, int n, const ms_real *range)
{
  if (i <= 0) {
    return 1.5*range[0]-0.5*range[1];
  }
  else if (i >= n) {
    return 1.5*range[n-1]-0.5*range[n-2];
  }
  else {
    return 0.5*(range[i-1]+range[i]);
  }
}


/* Simple single scattering calculation with no delta-Eddington
   scaling, usually used in conjunction with the TDTS multiple
   scattering algorithm */
int
ms_singlescatter(
    /* Input data */
    int n,                    /* Number of range gates in profile */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *ext,       /* Cloud/aerosol extinction coefficient, m-1 */
    const ms_real *ext_bscat_ratio,/* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft, m-1 (NULL for vacuum) */
    /* Output data */
    ms_real *bscat_out,       /* Measured backscatter, m-1 sr-1 */ 
    ms_real *bscat_air_out)
{
  ms_real optical_depth = 0.0;
  ms_real transmittance = 1.0;
  ms_real bscat_ext_ratio_air = 3.0 / (8.0 * MS_PI);
  int i;

  /* Assume that for wavelengths longer than 1 micron the gaseous
     extinction is entirely due to absorption rather than Rayleigh
     scattering */
  if (instrument.wavelength > 1e-6) {
    bscat_ext_ratio_air = 0.0;
  }

  if (ext_air) {
    /* Scattering atmospheric gases are present */
    for (i = 0; i < n; i++) {
      /* optical depth of layer */
      ms_real two_drange = 2.0*ms_get_drange(i, n, range);
      ms_real optical_depth_layer = (ext_air[i]+ext[i])*two_drange;
      ms_real factor;
      if (optical_depth_layer > 0.0) {
	ms_real transmittance_next;
	optical_depth += optical_depth_layer;
	transmittance_next = exp(-optical_depth);
	factor = (transmittance-transmittance_next)/optical_depth_layer;
	transmittance = transmittance_next;
      }
      else {
	factor = transmittance;
      }

      if (bscat_air_out) {
	/* Separate particulate and molecular returns */
	bscat_out[i] = factor*ext[i]/ext_bscat_ratio[i];
	bscat_air_out[i] = factor*ext_air[i]*bscat_ext_ratio_air;
      }
      else {
	bscat_out[i] = factor*(ext[i]/ext_bscat_ratio[i] 
			       + ext_air[i]*bscat_ext_ratio_air);
      }
    }
  }
  else {
    /* No atmospheric gases to consider for scattering, although they
       might attenuate */
    for (i = 0; i < n; i++) {
      ms_real two_drange = 2.0*ms_get_drange(i, n, range);
      ms_real optical_depth_layer;
      if (ext_air) {
	optical_depth_layer = (ext_air[i]+ext[i])*two_drange;
      }
      else {
	optical_depth_layer = ext[i]*two_drange;
      }

      if (optical_depth_layer > 0.0) {
	ms_real transmittance_next;
	optical_depth += optical_depth_layer;
	transmittance_next = exp(-optical_depth);
	bscat_out[i] = (transmittance-transmittance_next)*ext[i]
	  / (optical_depth_layer*ext_bscat_ratio[i]);
	transmittance = transmittance_next;
      }
      else {
	bscat_out[i] = 0.0;
      }
    }
    if (bscat_air_out) {
      for (i = 0; i < n; i++) {
	bscat_air_out[i] = 0.0;
      }
    }
  }

  return MS_SUCCESS;
}


