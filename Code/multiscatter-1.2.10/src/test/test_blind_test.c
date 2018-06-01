#include <gsl/gsl_multimin.h>
#include "multiscatter.h"

typedef struct {
  /* Input */
  ms_real* range;
  ms_real* radius;
  ms_real* ssa;
  ms_real* g;
  ms_real* ext_air;
  ms_real* ssa_air;
  ms_real* ext_bscat_ratio;
  ms_instrument* instrument;
  ms_surface* surface;
  ms_config* config;
  ms_real* bscat;
  ms_real* bscat_air;
  int n;
  int nx;
  int offset;
  /* Output */
  ms_real* bscat_fwd;
  ms_real* bscat_air_fwd;

  /* Adjoints */
  ms_real* bscat_AD;
  ms_real* bscat_air_AD;
  ms_real* ssa_AD;
  ms_real* g_AD;
  ms_real* ext_bscat_ratio_AD;

  /* Errors */
  ms_real* bscat_lninverrvar;
  ms_real* bscat_air_lninverrvar;

} ms_data;

#define EXT_FACTOR 10000.0
ms_real last_J = 0.0;
int num_J = 0;
int num_delJ = 0;

#define USE_LOG 1
//#define RETRIEVE_EXT_BSCAT_RATIO
//#define BSCAT_PHOTON 1e-6

ms_real
calc_J(const gsl_vector* v, void* params)
{
  ms_data* data = (ms_data*) params;
  int i, status;
  ms_real myJ = 0.0;
  ms_real ext[data->n];
#ifdef RETRIEVE_EXT_BSCAT_RATIO
  ms_real ext_bscat_ratio = exp(gsl_vector_get(v, data->nx));
#endif

  for (i = 0; i < data->n; i++) {
    ext[i] = 0.0;
  }

  for (i = 0; i < data->nx; i++) {
#ifdef USE_LOG
    ms_real ext_tmp = exp(gsl_vector_get(v, i));
#else
    ms_real ext_tmp = gsl_vector_get(v, i);
#endif
    if (ext_tmp <= 0.0 || isnan(ext_tmp)) {
      ext[i+data->offset] = 1.0e-6;
      myJ += ext_tmp*ext_tmp*EXT_FACTOR;
      printf(">>> #");
    }
    else {
      ext[i+data->offset] = ext_tmp;
    }
    fprintf(stderr, "?%14.8g ", ext_tmp);
#ifdef RETRIEVE_EXT_BSCAT_RATIO
    data->ext_bscat_ratio[i+data->offset] = ext_bscat_ratio;
#endif
  }

  /*  
  for (i = 0; i < data->n; i++) {
  fprintf(stderr, "range=%g radius=%g ext=%g ssa=%g g=%g ext_bscat_ratio=%g\n"
	  "ext_air=%g ssa_air=%g\n",
	  data->range[i], data->radius[i], ext[i], data->ssa[i],
	  data->g[i], data->ext_bscat_ratio[i], data->ext_air[i],
	  data->ssa_air[i]);
  }
  */
  status = multiscatter(data->n, data->n, data->config, 
			*(data->instrument), *(data->surface),
			data->range, data->radius, ext,
			data->ssa, data->g, data->ext_bscat_ratio,
			data->ext_air, data->ssa_air, NULL, NULL, 
			data->bscat_fwd, data->bscat_air_fwd);
  //  fprintf(stderr, "\n");
  for (i = 0; i < data->nx; i++) {
    ms_real diff = data->bscat_fwd[i+data->offset]
      - data->bscat[i+data->offset];
    /* Assume error variance (on the denominator) is proportional to
       backscatter (Poisson stats) */
    //    fprintf(stderr, "%g %g %g\n", v->data[i], data->bscat_fwd[i], data->bscat[i]);
    myJ += diff*diff * data->bscat_lninverrvar[i+data->offset];
    diff = data->bscat_air_fwd[i+data->offset]
      - data->bscat_air[i+data->offset];
    myJ += diff*diff * data->bscat_air_lninverrvar[i+data->offset];
    //    fprintf(stderr, " ---- %14.8g %14.8g %14.8g %14.8g\n",
    //	    data->bscat[i], data->bscat_fwd[i],
    //	    data->bscat_air[i], data->bscat_air_fwd[i]);
  }
  //  fprintf(stderr, "%g ", last_J - myJ);
  //    if (J == last_J) {
  //      J*=2;
  //    }

  if (last_J == myJ) {
    //    printf("%10.3g %10.3g %10.3g -> %g\n", ext[0], ext[6],
    //	   ext[12], J);
    fprintf(stderr,"*");
  }

  last_J = myJ;
  myJ *= 0.5;

  fprintf(stderr, "%0.18g\n", myJ);
  num_J++;
  return myJ;
}

void
calc_J_delJ(const gsl_vector* v, void* params, ms_real *myJ, gsl_vector* delJ)
{
  ms_data* data = (ms_data*) params;
  int status;
  int i;
  ms_real ext[data->n];
  ms_real ext_AD[data->n];
  ms_real ext_bscat_ratio_ADsum = 0.0;

  *myJ = calc_J(v, params); /* Needed to get bscat_fwd */

  for (i = 0; i < data->n; i++) {
    ext_AD[i] = data->ssa_AD[i] = data->g_AD[i]
      = data->ext_bscat_ratio_AD[i] = 0.0;
    data->bscat_AD[i] = data->bscat_air_AD[i] = 0.0;
    ext[i] = 0.0;
  }

  for (i = 0; i < data->nx; i++) {
#ifdef USE_LOG
    ms_real ext_tmp = exp(gsl_vector_get(v, i));
#else
    ms_real ext_tmp = gsl_vector_get(v, i);
#endif
    ext_AD[i+data->offset] = 0.0;
    if (ext_tmp <= 0.0) {
      ext[i+data->offset] = 1.0e-6;
      ext_AD[i+data->offset] += EXT_FACTOR*ext_tmp;
    }
    else {
      ext[i+data->offset] = ext_tmp;
    }
    data->bscat_AD[i+data->offset] 
      = (data->bscat_fwd[i+data->offset] - data->bscat[i+data->offset]) 
      * data->bscat_lninverrvar[i+data->offset];
    data->bscat_air_AD[i+data->offset] 
      = (data->bscat_air_fwd[i+data->offset] - data->bscat_air[i+data->offset])
      * data->bscat_lninverrvar[i+data->offset];
  }
  //  fprintf(stderr,"!");
  status = multiscatter_AD(data->n, data->n, data->config, 
			   *(data->instrument), *(data->surface),
			   data->range, data->radius, ext,
			   data->ssa, data->g, data->ext_bscat_ratio,
			   data->ext_air, data->ssa_air, NULL, NULL, 
			   data->bscat_fwd, data->bscat_air_fwd,
			   data->bscat_AD, data->bscat_air_AD,
			   ext_AD, data->ssa_AD, data->g_AD,
			   data->ext_bscat_ratio_AD);
  //  i=10;
  /*
  for (i = 0; i < data->n; i++) {
  fprintf(stderr, "range=%g radius=%g ext=%g ssa=%g g=%g ext_bscat_ratio=%g\n"
	  "ext_air=%g ssa_air=%g bscat_fwd=%g bscat_AD=%g\n",
	  data->range[i], data->radius[i], ext[i], data->ssa[i],
	  data->g[i], data->ext_bscat_ratio[i], data->ext_air[i],
	  data->ssa_air[i], data->bscat_fwd[i], data->bscat_AD[i]);
  }
  */
  
  for(i = -1; i <= data->n; i++) {
    fprintf(stderr, "%14.8g ", data->bscat_AD[i]);
  }
  fprintf(stderr, "bscat_AD\n");
  for(i = -1; i <= data->n; i++) {
    fprintf(stderr, "%14.8g ", data->bscat_air_AD[i]);
  }
  fprintf(stderr, "bscat_air_AD\n");

  for(i = 0; i < data->n; i++) {
    fprintf(stderr, "%14.8g ", data->ssa_AD[i]);
  }
  fprintf(stderr, "ssa_AD\n");

  for(i = 0; i < data->n; i++) {
    fprintf(stderr, "%14.8g ", data->g_AD[i]);
  }
  fprintf(stderr, "g_AD\n");
  for(i = 0; i < data->n; i++) {
    fprintf(stderr, "%14.8g ", data->ext_bscat_ratio_AD[i]);
  }
  fprintf(stderr, "ext_bscat_ratio_AD\n");
  for(i = 0; i < data->n; i++) {
    fprintf(stderr, "%14.8g ", ext[i]);
  }
  fprintf(stderr, "ext\n");
  for(i = 0; i < data->nx; i++) {
    fprintf(stderr, "%14.8g ", ext_AD[i+data->offset]);
    if (isnan(ext_AD[i+data->offset])) {
      ext_AD[i+data->offset] = 0.0;
    }
#ifdef USE_LOG
    gsl_vector_set(delJ, i, ext_AD[i+data->offset]*ext[i+data->offset]);
#else
    gsl_vector_set(delJ, i, ext_AD[i+data->offset]);
#endif
    //    fprintf(stderr,"ext_AD[%d]=%g ext=%g bscat_fwd=%g delJ=%g\n", i, ext_AD[i], ext[i],
    //    	    data->bscat_fwd[i], gsl_vector_get(delJ, i));
    ext_bscat_ratio_ADsum += data->ext_bscat_ratio_AD[i+data->offset]
      *data->ext_bscat_ratio[i+data->offset];
  }
#ifdef RETRIEVE_EXT_BSCAT_RATIO
  gsl_vector_set(delJ, data->nx, ext_bscat_ratio_ADsum);
#endif

  ext_AD[0] = 0.0;
  fprintf(stderr, "ext_AD\n");
  num_delJ++;
  return;
}

void
calc_delJ(const gsl_vector* v, void* params, gsl_vector* delJ)
{
  ms_real myJ;
  calc_J_delJ(v, params, &myJ, delJ);
  return;
}





int
main(int argc, char** argv)
{
  FILE* file;
  int n = 201;
  ms_real range[n];
  ms_real radius[n];
  ms_real ext[n];
  ms_real ssa[n];
  ms_real g[n];
  ms_real bscat[n];
  ms_real bscat_fwd[n];
  ms_real bscat_air[n];
  ms_real bscat_air_fwd[n];
  ms_real bscat_error[n];
  ms_real bscat_air_error[n];
  ms_real bscat_lninverrvar[n];
  ms_real bscat_air_lninverrvar[n];
  ms_real ext_air[n];
  ms_real ssa_air[n];
  ms_real ext_bscat_ratio[n];
  ms_real bscat_AD[n];
  ms_real bscat_air_AD[n];
  ms_real ssa_AD[n];
  ms_real g_AD[n];
  ms_real ext_bscat_ratio_AD[n];

  ms_instrument instrument = MS_DEFAULT_INSTRUMENT;
  ms_surface surface = MS_DEFAULT_SURFACE;
  ms_config config = MS_DEFAULT_CONFIG;
  ms_real rho_receiver[1] = { 0.0001 };
  ms_data data = {range, radius, ssa, g, ext_air, ssa_air, 
		  ext_bscat_ratio, &instrument, &surface,
		  &config, bscat, bscat_air, 
		  n, 0, 0, bscat_fwd, bscat_air_fwd, 
		  bscat_AD, bscat_air_AD, ssa_AD, g_AD,
		  ext_bscat_ratio_AD, bscat_lninverrvar,
		  bscat_air_lninverrvar};

  int i, status;

  const gsl_multimin_fdfminimizer_type *minimizer_type;
  gsl_multimin_fdfminimizer *minimizer;
  gsl_multimin_function_fdf function;
  gsl_vector* x;
  gsl_vector* delJ;
  size_t iter = 0;
  ms_real myJ = 0.0;
  ms_real myJnew = 0.0;

  minimizer_type = gsl_multimin_fdfminimizer_conjugate_pr;
  minimizer_type = gsl_multimin_fdfminimizer_vector_bfgs2;
  //      minimizer_type = gsl_multimin_fdfminimizer_steepest_descent;

  function.n = 0;
  function.f = calc_J;
  function.df = calc_delJ;
  function.fdf = calc_J_delJ;
  function.params = &data;

  data.range = range;
  data.radius = radius;
  data.ssa = ssa;
  data.g = g;

  file = fopen("BT_profile_11_60m_earthcare.dat", "r");


  for (i = 0; i < n; i++) {
    ms_real bscat_air_unatten;
    fscanf(file, "%lg %lg %lg %lg %lg %lg %lg %lg\n",
	   range+i, ext+i, ext_bscat_ratio+i,
	   &bscat_air_unatten,
	   bscat+i, bscat_air+i, bscat_error+i,
	   bscat_air_error+i);
    ext_air[i] = bscat_air_unatten*(8.0*M_PI/3.0);

    /*    
    range[i] = 10000.0-100.0*i;//*100.0;
    ext[i] = 1.0e-4;*/
    ssa[i] = 0.99;
    g[i] = 0.7;
    //    ext_air[i] = 1.0e-4;
    ssa_air[i] = 1.0;
    //    bscat[i] = bscat_air[i] = 0.0;
    //    ext_bscat_ratio[i] = 20.0;
    radius[i] = 1.0e-3;
    ssa_AD[i] = g_AD[i] = ext_bscat_ratio_AD[i] = 0.0;
    bscat_AD[i] = bscat_air_AD[i] = 0.0;


    if (bscat[i] > 0.0 && !data.offset) {
      data.offset = i;
    }
    else if (bscat[i] <= 0.0 && data.offset && !data.nx) {
      data.nx = i - data.offset;
    }

    bscat_lninverrvar[i] = 1.0
      /(bscat_error[i]*bscat_error[i]);
    bscat_air_lninverrvar[i] = 1.0e-2
      /(bscat_air_error[i]*bscat_air_error[i]);
    if (isinf(bscat_lninverrvar[i])) { 
      bscat_lninverrvar[i] = 0.0;
    }
    if (isinf(bscat_air_lninverrvar[i])) { 
      bscat_air_lninverrvar[i] = 0.0;
    }

    fprintf(stderr, "%d %g %g\n", i, data.bscat[i], data.bscat_lninverrvar[i]);

  }
#ifdef RETRIEVE_EXT_BSCAT_RATIO
  function.n = data.nx+1;
#else
  function.n = data.nx;
#endif

  instrument.receiver_type = TOP_HAT;
  instrument.altitude = 450.0e3;
  instrument.wavelength = 355.0e-9;
  instrument.rho_transmitter = 0.0001;
  instrument.rho_receiver = rho_receiver;
  instrument.nfov = 1;

  config.options |= MS_QUIET | MS_CHECK_FOR_NAN;
  //  config.options |= MS_SINGLE_ONLY;
  config.options |= MS_QSA_ONLY | MS_FAST_QSA;

  /*
  status = multiscatter(n, n, &config, 
			instrument, surface,
			range, radius, ext,
			ssa, g, ext_bscat_ratio,
			ext_air, ssa_air, NULL, NULL, 
			bscat, bscat_air);
  */

  x = gsl_vector_alloc(function.n);
#ifdef USE_LOG
  gsl_vector_set_all(x, log(2.0e-4));
  //  gsl_vector_set(x, 0, log(ext[0]));
  //  gsl_vector_set(x, 1, log(ext[1]));
  //  gsl_vector_set(x, 2, log(ext[2]));
  //  gsl_vector_set(x, 3, log(ext[3]));
#else
  gsl_vector_set_all(x, 2.0e-4);
  gsl_vector_set(x, 0, ext[0]);
#endif

#ifdef RETRIEVE_EXT_BSCAT_RATIO
  gsl_vector_set(x, data.nx, log(20.0));
#endif

  /*
  delJ = gsl_vector_alloc(n);
  calc_J_delJ(x, &data, &myJ, delJ);
  i = 1;
  fprintf(stderr, "J=%g, delJ=%g\n", myJ, gsl_vector_get(delJ, i));

  //  calc_delJ(x, &data, delJ);
  //  fprintf(stderr, "J=%g, delJ=%g\n", myJ, gsl_vector_get(delJ, i));

  gsl_vector_set(x, i, 2.2e-4);
  calc_J_delJ(x, &data, &myJnew, delJ);
  //  myJnew = calc_J(x, &data);
  fprintf(stderr, "J=%g, delJ=%g\n", myJnew, (myJnew-myJ)/0.2e-4);
  return 0;
  */

  //  for (i = 5; i < 15; i++) {
    //    gsl_vector_set(x, i, 1.0e-3);
  //  }
  calc_J(x, &data);

  minimizer = gsl_multimin_fdfminimizer_alloc(minimizer_type, function.n);
#ifdef USE_LOG
  gsl_multimin_fdfminimizer_set(minimizer, &function, x, 1.0, 0.9);
#else
  gsl_multimin_fdfminimizer_set(minimizer, &function, x, 0.0001, 0.1);
#endif


  fprintf(stderr, "*** FIRST GUESS ***\n");
  calc_J(x, &data);
  fprintf(stderr, "*******************\n");


  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(minimizer);
    if (status) {
      /*
      gsl_multimin_fdfminimizer_restart(minimizer);
      for (i = 0; i < n; i++) {
	if (gsl_vector_get(x, i) < 0.0) {
	  gsl_vector_set(x, i, 0.0);
	}
      }
      gsl_multimin_fdfminimizer_set(minimizer, &function, x, 0.0001, 0.1);
      */
      printf("@@@@@@@@@@@@@@@ An error occurred, code=%d\n", status);
      break;
    }

    status = gsl_multimin_test_gradient(minimizer->gradient, 1e-2);
    if (minimizer->f < 0.1/data.nx) {
      status = GSL_CONTINUE+1;
    }
    printf(">>> %d %g\n", iter, minimizer->f);
    //, gsl_vector_get(minimizer->x, data.nx));
    for (i = 0; i < data.nx; i++) {
      printf("%2d %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g %10.3g\n", i, 
	     ext[data.offset+i], 
#ifdef USE_LOG
	     exp(gsl_vector_get(minimizer->x, i)),
#else
	     gsl_vector_get(minimizer->x, i),
#endif
	     bscat[data.offset+i], bscat_fwd[data.offset+i],
	     bscat_air[data.offset+i], bscat_air_fwd[data.offset+i],
	     bscat_lninverrvar[data.offset+i], bscat_air_lninverrvar[data.offset+i]);
    }
  }  
  while (status == GSL_CONTINUE && iter < 100);
  //while (iter < 100);
  /*
  printf("Final pass (status=%d)\n", status);
  for (i = 0; i < n; i++) {
    printf("%2d %10.3g %10.3g %10.3g %10.3g\n", i, ext[i], 
#ifdef USE_LOG
	   exp(gsl_vector_get(minimizer->x, i)),
#else
	   gsl_vector_get(minimizer->x, i),
#endif
	   bscat[i], bscat_fwd[i]);
  }
  */
  printf(">>> Number of calculations of J: %d, delJ: %d\n",
	  num_J, num_delJ);

  return 0;
}

