/* fortran_interface.c -- Fortran interface to lidar multiple
   scattering algorithm

   Copyright (C) 2004-2011 Robin Hogan <r.j.hogan@reading.ac.uk> 

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

#include <stdlib.h>

#include "ms.h"

/* Number of contexts currently stored */
static int _n_context = 0;

/* List of dynamically allocated contexts */
static ms_context** _ms_context = NULL;

/* Return 1 if the context exists, 0 otherwise */
static
inline
int
context_exists(int context)
{
  if (context < _n_context && context >= 0 && _ms_context[context]) {
    return 1;
  }
  else {
    return 0;
  }
}

/* Create a new context and return its ID, but don't populate it  */
static
int
create_new_context()
{
  if (_n_context) {
    /* Search for a gap in the existing list of contexts */
    int i;
    ms_context** new_context_list;
    for (i = 0; i < _n_context; i++) {
      if (_ms_context[i] == NULL) {
	_ms_context[i] = malloc(sizeof(ms_context));
	if (!_ms_context[i]) {
	  fprintf(stderr, "Failed to allocate memory for multiscatter context\n");
	  return -1;
	}
	else {
	  /* New context successfully created */
	  return i;
	}
      }
    }
    /* No gap in existing context list */
    new_context_list = realloc(_ms_context, sizeof(ms_context*)*_n_context*2);
    if (!new_context_list) {
      fprintf(stderr, "Failed to allocate memory to extend list of multiscatter contexts\n");
      return -1;
    }
    _ms_context = new_context_list;
    for (i = _n_context; i < _n_context*2; i++) {
      _ms_context[i] = NULL;
    }
    i = _n_context;
    _n_context *= 2;
    _ms_context[i] = malloc(sizeof(ms_context));
    if (!_ms_context[i]) {
      fprintf(stderr, "Failed to allocate memory for multiscatter context\n");
      return -1;
    }
    else {
      /* New context successfully created */
      return i;
    }
  }

  _ms_context = malloc(sizeof(ms_context*));
  if (!_ms_context) {
    fprintf(stderr, "Failed to allocate memory to start list of multiscatter contexts\n");
    return -1;
  }
  _n_context = 1;
  _ms_context[0] = malloc(sizeof(ms_context));
  if (!_ms_context[0]) {
    fprintf(stderr, "Failed to allocate memory for multiscatter context\n");
    return -1;
  }
  else {
    /* New context successfully created */
    return 0;
  }
}

/* Fortran interface: Note that Fortran functions are stored in object
   files with a suffixed underscore. */

/* Free memory allocated to a particular context; this context can no
   longer be accessed */
void
ms_free_context_(int* context)
{
  if (context_exists(*context)) {
    _ms_context[*context] = NULL;
  }
}

/* Create a new context and return its integer ID */
int
ms_new_context_(ms_real* wavelength, ms_real* rho_transmitter,
		int* nfov, ms_real* rho_receiver)
{
  ms_context default_context = MS_DEFAULT_CONTEXT;

  ms_instrument* instrument;
  int context = create_new_context();
  int i;
  if (context < 0) {
    return context;
  }
  
  *(_ms_context[context]) = default_context;
  instrument = &(_ms_context[context]->instrument);
  instrument->wavelength = *wavelength;
  instrument->rho_transmitter = *rho_transmitter;
  instrument->nfov = *nfov;
  instrument->rho_receiver = malloc(sizeof(ms_real)* (*nfov));
  if (!instrument->rho_receiver) {
    fprintf(stderr, "Failed to allocate memory for receiver fields-of-view\n");
    ms_free_context_(&context);
    return -1;
  }
  for (i = 0; i < *nfov; i++) {
    instrument->rho_receiver[i] = rho_receiver[i];
  }
  /* Set other values */
  if (*wavelength > MS_RADAR_LIDAR_TRANSITION_WAVELENGTH) {
    /* Default settings for radar (lidar is the default) */
    instrument->receiver_type = MS_GAUSSIAN;
    _ms_context[context]->config.small_angle_algorithm
      = MS_SINGLE_AND_SMALL_ANGLE_NONE;
    _ms_context[context]->config.wide_angle_algorithm
      = MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE;
  }

  return context;
}

void
ms_set_altitude_(int* context, ms_real* altitude)
{
  ms_instrument* instrument;
  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in subroutine ms_set_altitude\n");
    return;
  }
  instrument = &(_ms_context[*context]->instrument);
  instrument->altitude = *altitude;
}


void
ms_set_gaussian_receiver_(int* context)
{
  ms_instrument* instrument;
  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in subroutine ms_set_gaussian_receiver\n");
    return;
  }
  instrument = &(_ms_context[*context]->instrument);
  instrument->receiver_type = MS_GAUSSIAN;
}

void
ms_set_top_hat_receiver_(int* context)
{
  ms_instrument* instrument;
  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in subroutine ms_set_top_hat_receiver\n");
    return;
  }
  instrument = &(_ms_context[*context]->instrument);
  instrument->receiver_type = MS_TOP_HAT;
}

void
ms_set_coherent_backscatter_enhancement_(int* context, ms_real* cbe)
{
  ms_config* config;
  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in subroutine ms_set_coherent_backscatter_enhancement\n");
    return;
  }
  config = &(_ms_context[*context])->config;
  config->coherent_backscatter_enhancement = *cbe;
}

/* Set the small-angle and wide-angle algorithms to be used according
   to the following:
    small_angle_code:
     0) No single-scattering or small-angle multiple scattering
     1) Single-scattering only
     2) Original small-angle algorithm (Hogan 2006)
     3) Fast small-angle algorithm (Hogan 2008); DEFAULT
     4) Explicit small-angle algorithm
    wide_angle_code:
     0) No wide-angle scattering
     1) Time-dependent two-stream algorithm (Hogan and Battaglia 2008): DEFAULT
*/
void
ms_set_algorithms_(int* context, int* small_angle_code, int* wide_angle_code)
{
  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in subroutine ms_set_small_angle_algorithm\n");
    return;
  }

  if (*small_angle_code < 0
      || *small_angle_code >= MS_NUM_SMALL_ANGLE_ALGORITHMS) {
    fprintf(stderr, "Invalid code in subroutine ms_set_small_angle_algorithm: must be 0 (no single scattering or small-angle multiple scattering), 1 (single scattering), 2 (Hogan 2006), 3 (faster Hogan 2008 algorithm) or 4 (explicit algorithm) but %d entered\n", 
	    *small_angle_code);
  }
  else {
    _ms_context[*context]->config.small_angle_algorithm = *small_angle_code;
  }

  if (*wide_angle_code < 0
      || *wide_angle_code >= MS_NUM_WIDE_ANGLE_ALGORITHMS) {
    fprintf(stderr, "Invalid wide-angle algorithm in subroutine ms_set_algorithms: must be 0 (no wide-angle scattering) or 1 (Hogan and Battaglia 2008) but %d entered\n", 
	    *wide_angle_code);
  }
  else {
    _ms_context[*context]->config.wide_angle_algorithm = *wide_angle_code;
  }
}

/* Select explicit small-angle scattering and choose the maximum
   number of scattering orders */
void
ms_set_explicit_max_scattering_order_(int* context, int* max_scattering_order)
{
  if (!context_exists(*context)) {
    fprintf(stderr, "Error: Invalid context in subroutine ms_set_explicit_max_scattering_order\n");
    return;
  }

  if (*max_scattering_order > 0) {
    _ms_context[*context]->config.max_scattering_order = *max_scattering_order;
    _ms_context[*context]->config.small_angle_algorithm
      = MS_SMALL_ANGLE_PVC_EXPLICIT;
  }
  else {
    fprintf(stderr, "Error: Max scattering order must be greater than 0 in subroutine ms_set_explicit_max_scattering_order but %d entered\n", *max_scattering_order);
  }
}

/* Multiple fields-of-view are treated as annular rather than overlapping */
void
ms_set_annular_detectors(int* context)
{
  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in subroutine ms_set_annular_detectors\n");
    return;
  }
  _ms_context[*context]->config.options |= MS_ANNULAR_DETECTORS;
}

void
ms_optimize_wide_angle_gates(int* context)
{
  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in subroutine ms_optimize_wide_angle_gates\n");
    return;
  }
  _ms_context[*context]->config.first_wide_angle_gate 
    = MS_AUTO_FIRST_WIDE_ANGLE_GATE;
}


/* Don't report anything except errors to standard error */
void
ms_set_quiet(int* context)
{
  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in subroutine ms_set_annular_detectors\n");
    return;
  }
  _ms_context[*context]->config.options |= MS_QUIET;
}

int
multiscatter_(
    /* Input data */
    int* context,
    int* n,                   /* Number of input gates */
    int* m,                   /* Number of output gates (>= n) */
    ms_real* range,           /* Height of each range gate, metres */
    ms_real* radius,          /* Particle equivalent radius, metres */
    ms_real* ext,             /* Particle extinction coefficient, m-1 */
    ms_real* ssa,             /* Particle single-scatter albedo */
    ms_real* g,               /* Particle asymmetry factor */
    ms_real* ext_bscat_ratio, /* Particle ext./backscatter ratio, sr */
    ms_real* ext_air,         /* Air ext. coefft., m-1 */
    ms_real* ssa_air,         /* Air single-scatter albedo */
    ms_real* droplet_fraction,     /* Fraction of extinction from droplets */
    ms_real* pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Output data */
    ms_real* bscat_out)       /* Measured backscatter, m-1 sr-1 */
{
  ms_config* config = &(_ms_context[*context]->config);
  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in function multiscatter\n");
    return MS_INVALID_CONTEXT;
  }
  
  if (!(config->options & MS_QUIET)) {
    ms_print_algorithms(*config, _ms_context[*context]->instrument,
			range, 0, 0, 0, 0, 0, stderr);
  }

  int status = multiscatter(*n, *m, config,
			    _ms_context[*context]->instrument,
			    _ms_context[*context]->surface,
			    range, radius, ext, ssa, g,
			    ext_bscat_ratio, ext_air, ssa_air,
			    droplet_fraction, pristine_ice_fraction,
			    bscat_out, NULL);
  if (status != MS_SUCCESS) {
    fprintf(stderr, "Error occurred in function multiscatter\n");
  }
  return status;
}


int
multiscatter_hsrl_(
    /* Input data */
    int* context,
    int* n,                   /* Number of input gates */
    int* m,                   /* Number of output gates (>= n) */
    ms_real* range,           /* Height of each range gate, metres */
    ms_real* radius,          /* Cloud/aerosol equivalent radius, microns */
    ms_real* ext,             /* Total extinction coefficient, m-1 */
    ms_real* ssa,             /* Total single-scatter albedo */
    ms_real* g,               /* Total asymmetry factor */
    ms_real* ext_bscat_ratio, /* Cloud/aerosol ext./backscatter ratio, sr */
    ms_real* ext_air,         /* Air ext. coefft., m-1 */
    ms_real* ssa_air,         /* Air single-scatter albedo */
    ms_real* droplet_fraction,     /* Fraction of extinction from droplets */
    ms_real* pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Output data */
    ms_real* bscat_out,       /* Measured backscatter, m-1 sr-1 */
    ms_real* bscat_air_out)   /* Measured backscatter of air, m-1 sr-1 */
{
  ms_config* config = &(_ms_context[*context]->config);
  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in function multiscatter_hsrl\n");
    return MS_INVALID_CONTEXT;
  }
  
  if (!(config->options & MS_QUIET)) {
    ms_print_algorithms(*config, _ms_context[*context]->instrument,
			range, 0, 1, 0, 0, 0, stderr);
  }

  int status = multiscatter(*n, *m, config,
			    _ms_context[*context]->instrument,
			    _ms_context[*context]->surface,
			    range, radius, ext, ssa, g,
			    ext_bscat_ratio, ext_air, ssa_air,
			    droplet_fraction, pristine_ice_fraction,
			    bscat_out, bscat_air_out);
  if (status != MS_SUCCESS) {
    fprintf(stderr, "Error occurred in function multiscatter_hsrl\n");
  }
  return status;
}

/* Just call small-angle part only thereby not needing
   single-scattering albedo (ssa), asymmetry factor (g), and we assume
   that this is a short wavelength lidar so don't need ssa_air (assume
   it is 1) */
int
multiscatter_small_angle_only_(
    /* Input data */
    int* context,
    int* n,                   /* Number of input gates */
    ms_real* range,           /* Height of each range gate, metres */
    ms_real* radius,          /* Cloud/aerosol equivalent radius, microns */
    ms_real* ext,             /* Total extinction coefficient, m-1 */
    ms_real* ext_bscat_ratio, /* Cloud/aerosol ext./backscatter ratio, sr */
    ms_real* ext_air,         /* Air ext. coefft., m-1 */
    ms_real* droplet_fraction,     /* Fraction of extinction from droplets */
    ms_real* pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Output data */
    ms_real* bscat_out)       /* Measured backscatter, m-1 sr-1 */
{
  ms_config* config = &(_ms_context[*context]->config);

  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in function multiscatter_small_angle_only\n");
    return MS_INVALID_CONTEXT;
  }
  
  if (!(config->options & MS_QUIET)) {
    ms_print_algorithms(*config, _ms_context[*context]->instrument,
			range, 0, 0, 0, 0, 0, stderr);
  }

  int status = multiscatter(*n, *n, config,
			    _ms_context[*context]->instrument,
			    _ms_context[*context]->surface,
			    range, radius, ext, NULL, NULL,
			    ext_bscat_ratio, ext_air, NULL,
			    droplet_fraction, pristine_ice_fraction,
			    bscat_out, NULL);
  if (status != MS_SUCCESS) {
    fprintf(stderr, "Error occurred in function multiscatter_small_angle_only\n");
  }
  return status;
}


/* Just call small-angle part only thereby not needing
   single-scattering albedo (ssa), asymmetry factor (g), and we assume
   that this is a short wavelength lidar so don't need ssa_air (assume
   it is 1) */
int
multiscatter_small_angle_only_hsrl_(
    /* Input data */
    int* context,
    int* n,                   /* Number of input gates */
    ms_real* range,           /* Height of each range gate, metres */
    ms_real* radius,          /* Cloud/aerosol equivalent radius, microns */
    ms_real* ext,             /* Total extinction coefficient, m-1 */
    ms_real* ext_bscat_ratio, /* Cloud/aerosol ext./backscatter ratio, sr */
    ms_real* ext_air,         /* Air ext. coefft., m-1 */
    ms_real* droplet_fraction,     /* Fraction of extinction from droplets */
    ms_real* pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Output data */
    ms_real* bscat_out,       /* Measured backscatter, m-1 sr-1 */
    ms_real* bscat_air_out)   /* Measured backscatter of air, m-1 sr-1 */
{
  ms_config* config = &(_ms_context[*context]->config);

  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in function multiscatter_small_angle_only_hsrl\n");
    return MS_INVALID_CONTEXT;
  }
  
  if (!(config->options & MS_QUIET)) {
    ms_print_algorithms(*config, _ms_context[*context]->instrument,
			range, 0, 1, 0, 0, 0, stderr);
  }

  int status = multiscatter(*n, *n, config,
			    _ms_context[*context]->instrument,
			    _ms_context[*context]->surface,
			    range, radius, ext, NULL, NULL,
			    ext_bscat_ratio, ext_air, NULL,
			    droplet_fraction, pristine_ice_fraction,
			    bscat_out, bscat_air_out);
  if (status != MS_SUCCESS) {
    fprintf(stderr, "Error occurred in function multiscatter_small_angle_only_hsrl\n");
  }
  return status;
}


/* Simple interface to multiscatter algorithm: small-angle multiple
   scattering only, phase function near backscatter is assumed
   isotropic so not droplet_fraction and pristine_ice_fraction, and
   finally the air and molecular backscatters are merged, so no
   bscat_air_out is required. */
int
multiscatter_simple_(
    /* Input data */
    int* context,
    int* n,                   /* Number of input gates */
    ms_real* range,           /* Height of each range gate, metres */
    ms_real* radius,          /* Cloud/aerosol equivalent radius, microns */
    ms_real* ext,             /* Total extinction coefficient, m-1 */
    ms_real* ext_bscat_ratio, /* Cloud/aerosol ext./backscatter ratio, sr */
    ms_real* ext_air,         /* Air ext. coefft., m-1 */
    /* Output data */
    ms_real* bscat_out)       /* Measured backscatter, m-1 sr-1 */
{
  ms_config* config = &(_ms_context[*context]->config);

  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in function multiscatter_simple\n");
    return MS_INVALID_CONTEXT;
  }
  
  if (!(config->options & MS_QUIET)) {
    ms_print_algorithms(*config, _ms_context[*context]->instrument,
			range, 1, 0, 0, 0, 0, stderr);
  }

  int status = multiscatter(*n, *n, config,
			    _ms_context[*context]->instrument,
			    _ms_context[*context]->surface,
			    range, radius, ext, NULL, NULL,
			    ext_bscat_ratio, ext_air, NULL,
			    NULL, NULL,
			    bscat_out, NULL);
  if (status != MS_SUCCESS) {
    fprintf(stderr, "Error occurred in function multiscatter_simple\n");
  }
  return status;
}


int
multiscatter_AD_(
    /* Input data */
    int* context,
    int* n,                   /* Number of input gates */
    int* m,                   /* Number of output gates (>= n) */
    ms_real* range,           /* Height of each range gate, metres */
    ms_real* radius,          /* Particle equivalent radius, metres */
    ms_real* ext,             /* Particle extinction coefficient, m-1 */
    ms_real* ssa,             /* Particle single-scatter albedo */
    ms_real* g,               /* Particle asymmetry factor */
    ms_real* ext_bscat_ratio, /* Particle ext./backscatter ratio, sr */
    ms_real* ext_air,         /* Air ext. coefft., m-1 */
    ms_real* ssa_air,         /* Air single-scatter albedo */
    ms_real* droplet_fraction,     /* Fraction of extinction from droplets */
    ms_real* pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Output data */
    ms_real* bscat_out,       /* Measured backscatter, m-1 sr-1 */
    /* Adjoint inputs */
    const ms_real *bscat_AD,
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
  ms_config* config = &(_ms_context[*context]->config);
  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in function multiscatter_AD\n");
    return MS_INVALID_CONTEXT;
  }
  
  if (!(config->options & MS_QUIET)) {
    ms_print_algorithms(*config, _ms_context[*context]->instrument,
			range, 0, 0, 0, 0, 0, stderr);
  }

  int status = multiscatter_AD(*n, *m, config,
			       _ms_context[*context]->instrument,
			       _ms_context[*context]->surface,
			       range, radius, ext, ssa, g,
			       ext_bscat_ratio, ext_air, ssa_air,
			       droplet_fraction, pristine_ice_fraction,
			       bscat_out, NULL,
			       bscat_AD, NULL,
			       radius_AD, ext_AD, ssa_AD, g_AD,
			       ext_bscat_ratio_AD, ext_air_AD,
			       droplet_fraction_AD, pristine_ice_fraction_AD);
  if (status != MS_SUCCESS) {
    fprintf(stderr, "Error occurred in function multiscatter_AD\n");
  }
  return status;
}

int
multiscatter_hsrl_AD_(
    /* Input data */
    int* context,
    int* n,                   /* Number of input gates */
    int* m,                   /* Number of output gates (>= n) */
    ms_real* range,           /* Height of each range gate, metres */
    ms_real* radius,          /* Particle equivalent radius, metres */
    ms_real* ext,             /* Particle extinction coefficient, m-1 */
    ms_real* ssa,             /* Particle single-scatter albedo */
    ms_real* g,               /* Particle asymmetry factor */
    ms_real* ext_bscat_ratio, /* Particle ext./backscatter ratio, sr */
    ms_real* ext_air,         /* Air ext. coefft., m-1 */
    ms_real* ssa_air,         /* Air single-scatter albedo */
    ms_real* droplet_fraction,     /* Fraction of extinction from droplets */
    ms_real* pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Output data */
    ms_real* bscat_out,       /* Measured particle backscatter, m-1 sr-1 */
    ms_real* bscat_air_out,   /* Measured air backscatter, m-1 sr-1 */
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
  ms_config* config = &(_ms_context[*context]->config);
  if (!context_exists(*context)) {
    fprintf(stderr, "Invalid context in function multiscatter_AD\n");
    return MS_INVALID_CONTEXT;
  }
  
  if (!(config->options & MS_QUIET)) {
    ms_print_algorithms(*config, _ms_context[*context]->instrument,
			range, 0, 0, 0, 0, 0, stderr);
  }

  int status = multiscatter_AD(*n, *m, config,
			       _ms_context[*context]->instrument,
			       _ms_context[*context]->surface,
			       range, radius, ext, ssa, g,
			       ext_bscat_ratio, ext_air, ssa_air,
			       droplet_fraction, pristine_ice_fraction,
			       bscat_out, bscat_air_out,
			       bscat_AD, bscat_air_AD,
			       radius_AD, ext_AD, ssa_AD, g_AD,
			       ext_bscat_ratio_AD, ext_air_AD,
			       droplet_fraction_AD, pristine_ice_fraction_AD);
  if (status != MS_SUCCESS) {
    fprintf(stderr, "Error occurred in function multiscatter_AD\n");
  }
  return status;
}
