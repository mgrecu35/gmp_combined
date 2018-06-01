/* multiscatter.c -- Quickbeam interface to TDTS multiple scattering method

   Copyright (C) 2008 Robin Hogan <r.j.hogan@reading.ac.uk> 

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


   The Quickbeam interface is provided by the multiscatter_()
   function; see the comments just before this function for details of
   how to call it. For further details about this algorithm, see the
   papers and other material here:
     http://www.met.rdg.ac.uk/clouds/multiscatter/
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

/* Note that this header file contains everything for the multiscatter
   distribution, only a subset of which is released with Quickbeam. */
#include "multiscatter.h"

/* Parameters that control behaviour of the algorithm */

/* For speed it is useful not to calculate multiple scattering on the
   profile, but only the part that is optically thick enough that
   multiple scattering will be significant. The threshold is in terms
   of the ratio of a mean-free-path to the radar footprint - if it is
   less than this value then multiple scattering will be calculated
   for all gates beyond this point. A typical number would be 20, but
   if it is set very large then multiple scattering will be calculated
   for the entire profile. */
#define MS_MFP_FOOTPRINT_THRESHOLD 20.0

/* The input data may be on an irregular grid, while the multiple
   scattering calculation is performed on a regular grid. In deciding
   the resolution to use for this regular grid, there is a trade-off
   between the accuracy and the speed. However, there is no need to
   run the algorithm at a higher resolution than the important initial
   region of the cloud, and after many optical depths into the cloud
   all sharp gradients are smoothed out anyway. To decide what
   resolution to use, the cloud is searched up to a particular
   scattering optical depth, and the finest grid spacing in this
   region is used for the regular grid. This variable determines the
   optical depth threshold to use. A typical value would be 5. */
#define MS_OD_DHEIGHT_THRESHOLD 5.0

/* Michael Mishchenko argues that coherent backscatter should enhance
   the multiple scattering by a factor of 2 for radar, but not for
   lidar */
#define MS_COHERENT_BACKSCATTER_ENHANCEMENT 2.0

/* To minimise the overhead associated with memory allocation, the
   memory allocated for performing the algorithm is kept in static
   variables so is available to the next call to the algorithm. If
   more space is required, it is resized upwards. Therefore, the
   algorithm will be faster per profile for larger numbers of
   profiles.  HOWEVER, THIS CODE IS NOT THREAD-SAFE. If
   multiscatter_() is to be called simultaneously by different threads
   of the same program then the memory allocation parts will need to
   be recoded */

/* Static variables only used in this file are prefixed with "msf_" to
   denote that they are associated with the multiple-scattering
   Fortran interface. These quantities are on the regular grid. */
static ms_real *msf_range = NULL;        /* Range (m) */
static ms_real *msf_ext = NULL;          /* Extinction coefficient (m-1) */
static ms_real *msf_ssa = NULL;          /* Single scattering albedo */
static ms_real *msf_g = NULL;            /* Asymmetry factor */
static ms_real *msf_src_power_in = NULL; /* Normalized source inwards power */
static ms_real *msf_src_power_out = NULL;/* Normalized source outwards power */
static ms_real *msf_src_width2 = NULL;   /* Squared width of source (m2) */
static ms_real *msf_bscat_out = NULL;    /* Apparent backscatter (m-1 sr-1) */

/* Number of gates currently allocated in the intermediate
   variables: */
static int msf_num_gates_allocated = 0;

/* Number of gates defined in the last call to the algorithm: */
static int msf_num_gates_defined = 0;

/* Allocate arrays used in this file, or use existing memory
   allocation if possible */
static
int
msf_init_intermediate_arrays(int num_gates)
{
  /* Macro for allocating or reallocating memory for a single array */
#define ALLOCATE(target, num_bytes) \
   if (target) { target--; } \
   if ((target = realloc(target, num_bytes)) == NULL) { target = NULL; \
     fprintf(stderr, \
       "%s line %d: Out of memory in attempt to allocate %d bytes\n", \
       __FILE__, __LINE__, num_bytes); exit(1); } else { target++; }


  int i;
  if (msf_num_gates_allocated < num_gates) {
    /* Number of gates allocated is less than number requested (which
       may be zero if this function has not been called): (re)allocate
       memory. */
    int num_bytes = (num_gates+2) * sizeof(ms_real);
    msf_num_gates_allocated = 0; /* In case error occurs */
    ALLOCATE(msf_range, num_bytes);
    ALLOCATE(msf_ext, num_bytes);
    ALLOCATE(msf_ssa, num_bytes);
    ALLOCATE(msf_g, num_bytes);
    ALLOCATE(msf_src_power_in, num_bytes);
    ALLOCATE(msf_src_power_out, num_bytes);
    ALLOCATE(msf_src_width2, num_bytes);
    ALLOCATE(msf_bscat_out, num_bytes);
  }
  msf_num_gates_allocated = num_gates;
  /* Set each element of each array to zero. */
  for (i = -1; i <= num_gates; i++) {
    msf_range[i] = 0.0;
    msf_ext[i] = 0.0;
    msf_ssa[i] = 0.0;
    msf_g[i] = 0.0;
    msf_src_power_in[i] = 0.0;
    msf_src_power_out[i] = 0.0;
    msf_src_width2[i] = 1.0;
    msf_bscat_out[i] = 0.0;
  }
  msf_num_gates_defined = num_gates;
  return MS_SUCCESS;
}

/* Options are held in this variable as a bit mask - see
   multiscatter.h for possibilities (note that many are for the
   small-angle part of the multiscatter distribution so will have no
   effect on the wide angle part used here */
int ms_options = 0;

/* Interface to wide-angle multiple scattering code for the Quickbeam
   radar simulator. This should be called from Fortran as follows:
      call multiscatter(...)
   Note that they are passed into C as pointers.

   This function interpolates the input data to a regular grid, it
   calls the time-dependent two-stream (TDTS) algorithm (the function
   multiscatter_tdts), and it interpolates the results back on to the
   input grid. */
void
multiscatter_(/* Scalar inputs */
	      ms_real* freq,        /* Radar frequency (GHz) */
	      ms_real* altitude,    /* Instrument altitude (km) */
	      ms_real* beamwidth,   /* 1/e full-angle one-way beamwidth (radians) */
	      ms_real* K2_ref,      /* Reference |K^2| (e.g. 0.75) */
	      int* ngate,           /* Number of gates in the vectors */
	      /* Vector inputs */
	      ms_real* height,      /* Height of each gate (km) */
	      ms_real* ext,         /* Extinction coefficient (m-1) */
	      ms_real* ssa,         /* Single scattering albedo */
	      ms_real* g,           /* Asymmetry factor */
	      /* Output */
	      ms_real* Z_multiscat  /* Multi-scat. radar reflectivity (mm6 m-3) */
	      )
{
  /* Structures to pass to the TDTS algorithm */
  ms_instrument instrument;
  ms_surface surface;

  ms_real rho = *beamwidth * 0.5;  /* The 1/e half-angle one-way beamwidth */
  ms_real scat_od = 0.0;           /* Scattering optical depth */
  ms_real od = 0.0;                /* Optical depth */

  /* Heights/indices at grid boundaries used for interpolation */
  ms_real start_edge_height = 0.0, end_edge_height = 0.0;
  ms_real edge_height1 = 0.0, edge_height2 = 0.0;
  ms_real j1, j2;

  ms_real dheight = *altitude; /* Resolution of regular grid (km) */
  ms_real direction = +1.0;    /* +1 for ground-based, -1 for spaceborne */

  /* Scaling factors */
  ms_real ss_multiplier = MS_COHERENT_BACKSCATTER_ENHANCEMENT; 
  ms_real Z_factor;

  int n = *ngate;    /* Number of irregular gates */
  int m = 0;         /* Number of regular gates */
  int i;             /* Index to irregular grid */
  int istart;        /* First gate with significant multiple scattering */
  int j;             /* Index to regular grid */

  int status;        /* Output from function call */
  
  /* Set direction */
  if (height[1] < height[0]) {
    direction = -1.0;
  }

  /* Find start gate where multiple scattering likely to be significant */
  istart = n;
  for (i = 0; i < n; i++) {
    Z_multiscat[i] = 0.0; /* Set gate to zero */
    //    fprintf(stderr, "%i,%g,%g,%g \t", i, ext[i],ssa[i],g[i]);

    if (ext[i] > 0.0 && ssa[i] > 0.0) {
      /* The footprint radius at the altitude of the cloud */
      ms_real footprint = fabs(*altitude-height[i])*rho*1000.0;
      /* The relevant mean-free-path */
      ms_real mfp = 1.0/(ext[i]*ssa[i]*(1.0 - g[i]));
      if (mfp < footprint*MS_MFP_FOOTPRINT_THRESHOLD) {
	/* Threshold met: start multiple scattering calculation from
	   this gate */
	istart = i;
	break;
      }
    }
    /* Estimate optical depth to the start of the multiple scattering
       region: this will be used in calculating the source terms for
       the TDTS calculation later on */
    if (i == 0) {
      od += ext[i]*fabs(height[1] - height[0])*1000.0;
    }
    else if (i < n-1) {
      od += ext[i]*0.5*fabs(height[i+1] - height[i-1])*1000.0;
    }
  }

  if (istart >= n-1) {
    /* Cloud is not optically thick enough for significant multiple
       scattering (or only the lowest gate exceeds the threshold):
       return with all the elements of Z_multiscat set to zero */
    Z_multiscat[n-1] = 0.0;
    return;
  }

  /* Calculate the height of the edge of the first multiple scattering
     gate, treating the grid edges as half way between the elements of
     the height vector */ 
  if (istart == 0) {
    /* Extrapolate */
    start_edge_height = 1.5*height[0] - 0.5*height[1];
  }
  else {
    /* Interpolate */
    start_edge_height = 0.5*(height[istart-1] + height[istart]);
  }

  /* Find the smallest grid spacing in the first
     MS_OD_DHEIGHT_THRESHOLD scattering optical depths of the cloud
     where the mean-free-path criterion above is also satisfied: this
     grid spacing will become "dheight" */
  /* Far edge of first gate */
  edge_height2 = 0.5*(height[istart] + height[istart+1]);

  /* Initial value for dheight is the spacing for the first gate */
  dheight = fabs(start_edge_height - edge_height2);

  /* Initial scattering optical depth */
  scat_od = 1000*dheight*ext[istart]*ssa[istart];

  /* Loop through subsequent gates */
  for (i = istart+1; i < n; i++) {
    ms_real dheight_i;            /* Grid spacing of gate i */
    edge_height1 = edge_height2;  /* Near edge of gate i is far edge of gate i-1 */
    /* Calculate location of far edge of gate i */
    if (i < n-1) {
      /* Interpolate */
      edge_height2 = 0.5*(height[i] + height[i+1]);
    }
    else {
      /* Extrapolate */
      edge_height2 = 1.5*height[i] - 0.5*height[i-1];
    }
    /* Set the grid spacing of gate i */
    dheight_i = fabs(edge_height1-edge_height2);

    /* Is this smaller than dheight, and is there multiple-scattering
       cloud present? */
    if (dheight > dheight_i && ext[i] > 0.0 && ssa[i] > 0.0) {
      ms_real footprint = fabs(*altitude-height[i])*rho*1000.0;
      ms_real mfp = 1.0/(ext[i]*ssa[i]*(1.0 - g[i]));
      if (mfp < footprint*MS_MFP_FOOTPRINT_THRESHOLD) {
	/* Threshold met: set this to the new grid spacing */
	dheight = dheight_i;
      }
    }
    /* Increment the scattering optical depth */
    scat_od += 1000*dheight*ext[i]*ssa[i];
    if (scat_od > MS_OD_DHEIGHT_THRESHOLD) {
      /* We have reached our threshold scattering optical depth: stop
	 looking for more finely spaced gates */
      break;
    }
  }
  /* Calculate location of far edge of the final gate */
  end_edge_height = 1.5*height[n-1] - 0.5*height[n-2];

  /* Assign instrument properties */
  instrument.receiver_type = GAUSSIAN;
  instrument.altitude = *altitude * 1000.0;
  instrument.wavelength = MS_C/(*freq*1.0e9); /* m */
  instrument.rho_transmitter = rho;
  instrument.rho_receiver = rho;

  /* Multiplier to ensure that the result scales to the
     single-scattering return, accounting for transmiter/receiver
     overlap */
  ss_multiplier = MS_COHERENT_BACKSCATTER_ENHANCEMENT
    + instrument.rho_transmitter*instrument.rho_transmitter
    /(instrument.rho_receiver*instrument.rho_receiver);

  /* Scaling factor to convert apparent backscatter coefficient (m-1
     sr-1) into radar reflectivity factor (mm6 m-3) */
  Z_factor = 1.0e18 * pow(instrument.wavelength/MS_PI, 4.0) * 4.0 / *K2_ref;

  /* Calculate how many points will be required */
  m = ((int) (fabs(start_edge_height - end_edge_height)/dheight));

  /* Allocate intermediate arrays, or re-zero existing memory if it is
     already allocated from a previous call */
  msf_init_intermediate_arrays(m);

  /* Interpolate scattering properties on to a regular grid */

  /* Select the first gate of the irregular grid that will be used in
     the interpolation */
  if (istart > 0) {
    i = istart-1;
  }
  else {
    i = 0;
  }

  /* Loop through the elements of the regular grid */
  for (j = 0; j < m; j++) {
    /* Set the range in metres */
    msf_range[j] = 1000.0*(start_edge_height + direction*dheight*(0.5+j));
    /* Find first irregular gate i beyond regular gate j */
    if (direction > 0.0) {
      while (i < n && height[i]*1000.0 < msf_range[j]) {
	i++;
      }
    }
    else {
      while (i < n && height[i]*1000.0 > msf_range[j]) {
	i++;
      }
    }

    /* 1/e squared half width of outgoing beam */
    msf_src_width2[j] = rho * rho
      * (msf_range[j] - *altitude*1000.0)
      * (msf_range[j] - *altitude*1000.0);

    /* Calculate values on regular grid */
    if (i == 0) {
      /* Rather than extrapolate, just use the values from the first
	 gate of the irregular grid */
      msf_ext[j] = ext[i];
      msf_ssa[j] = ssa[i];
      msf_g[j] = g[i];
    }
    else {
      /* Perform linear interpolation using the following weights */
      ms_real weight1 
	= (height[i]-0.001*msf_range[j])/(height[i]-height[i-1]);
      ms_real weight2 
	= (0.001*msf_range[j]-height[i-1])/(height[i]-height[i-1]);
      msf_ext[j] = weight1*ext[i-1] + weight2*ext[i];

      /* Interpolate the scattering coefficient and the
	 asymmetry-weighted scattering coefficient */
      ms_real scat_coefft
	= weight1*ext[i-1]*ssa[i-1] + weight2*ext[i]*ssa[i];
      ms_real g_scat_coefft
	= weight1*ext[i-1]*ssa[i-1]*g[i-1] + weight2*ext[i]*ssa[i]*g[i];

      /* Calculate scattering properties on regular grid */
      msf_ssa[j] = scat_coefft / msf_ext[j];
      if (scat_coefft > 0.0) {
	msf_g[j] = g_scat_coefft / scat_coefft;
      }
    }

    /* If there is any scattering at the gate then calculate the
       source powers for the TDTS calculation*/
    if (msf_ext[j] > 0.0 && msf_ssa[j] > 0.0) {
      ms_real layer_od = msf_ext[j]*dheight*1000.0; /* Optical depth of gate */
      ms_real src_power = ss_multiplier*exp(-od)
	*(1-exp(-msf_ext[j]*msf_ssa[j]*dheight*1000.0));
      /* Use two-stream phase function to predicts how much of
	 initially scattered energy is into the forward (out) and
	 backward (in) hemisphere, but account for the fact that it
	 predicts negative values for asymmetry factor outside the
	 range -2/3 to 2/3 */
      if (msf_g[j] > 1.0/(3.0*MS_MU1)) {
	msf_src_power_in[j] = 0.0;
	msf_src_power_out[j] = src_power;
      }
      else if (msf_g[j] < -1.0/(3.0*MS_MU1)) {
	msf_src_power_in[j] = src_power;
	msf_src_power_out[j] = 0.0;
      }
      else {
	msf_src_power_in[j] = src_power*0.5*(1.0-3.0*msf_g[j]*MS_MU1);
	msf_src_power_out[j] = src_power*0.5*(1.0+3.0*msf_g[j]*MS_MU1);
      }
      od += layer_od; /* Increment the optical depth */
    }
  }

  /*
  for (i = 0; i < n; i++) {
    fprintf(stderr, "%3d %9.4g %9.4g %9.4g %9.4g\n",
	    i, height[i], ext[i], ssa[i], g[i]);
  }
  */

  fprintf(stderr, "Calculating radar multiple scattering for profile with optical depth %g\n",
	  od);

  /* Call the TDTS algorithm */
  status = multiscatter_tdts(m, m, instrument, surface, msf_range,
			     msf_ext, msf_ssa, msf_g,
			     msf_src_power_in, msf_src_power_out,
			     msf_src_width2, msf_bscat_out);
  /* Exit if an error occurred */
  if (status != MS_SUCCESS) {
    fprintf(stderr,
	    "%s line %d: error code %d returned from multiscatter_tdts function\n",
	    __FILE__, __LINE__, status);
    exit(status);
  }
  /*
  fprintf(stderr, "  i    height     ext      ssa         g      src-in   src-out  src-width2  bscat\n");
  for (j = 0; j < m; j++) {
    fprintf(stderr, "%3d %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g\n",
	    j, msf_range[j]*0.001, msf_ext[j], msf_ssa[j], msf_g[j],
	    msf_src_power_out[j], msf_src_power_in[j], msf_src_width2[j],
	    msf_bscat_out[j]);
  }
  */

  /* Now interpolate the backscatter back on to the irregular grid,
     but conserving the integral */

  /* Define j1 and j2 as the corresponding indices of the regular
     grid, but noting that they may not be integers; note that here we
     treat msf_bscat_out[0] as being a top-hat function between j=0
     and j=1. By design the first edge was at j=0*/
  j1 = 0;

  /* Loop through each gate of the irregular grid */
  for (i = istart; i < n; i++) {
    /* Set apparent backscatter coefficient of irregular gate i
       initially to zero */
    ms_real bscat = 0.0;
    int j;

    /* Calculate the far boundary, edge_height2 */
    if (i < n-1) {
      /* Interpolate */
      edge_height2 = 0.5*(height[i] + height[i+1]);
    }
    else {
      /* Extrapolate */
      edge_height2 = 1.5*height[i] - 0.5*height[i-1];
    }

    /* Calculate the corresponding equivalent index in j-space */
    j2 = fabs(start_edge_height - edge_height2)/dheight;

    /* Loop through the regular gates that contribute to the irregular
       gate; the first one to consider is the rounded-down j1 */
    j = (int) j1;
    /*    fprintf(stderr, "Gate %d: %g-%g\n", i, j1, j2); */
    while (j < j2) {
      /* Calculate how much of the regular gate contributes to the
	 irregular gate, by defining j_upper and j_lower, which
	 delimit the part of the j to j+1 range that contributes to
	 the irregular grid */
      ms_real j_upper = j+1.0;
      ms_real j_lower = j;
      if (j_upper > j2) {
	j_upper = j2;
      }
      if (j_lower < j1) {
	j_lower = j1;
      }
      /* fprintf(stderr, "   Alloc %d %g\n", j, (j_upper-j_lower)/(j2-j1)); */
      /* Increment backscatter at irregular gate i by the contribution
	 from regular gate j */
      bscat += msf_bscat_out[j] * (j_upper-j_lower)/(j2-j1);
      j++;
    }
    /* Convert apparent backscatter coefficient to effective radar
       reflectivity factor (due to multiple scattering only) */
    Z_multiscat[i] = bscat * Z_factor;

    /* Set the near boundary of the next gate to be the far boundary
       of the previous gate */
    j1 = j2;
  }

  return;
}
