#include "multiscatter.h"

typedef struct {
  ms_real power;
  ms_real power_x_width2;
  ms_real power_x_zeta2;
  ms_real power_x_cov;
} pvc_distribution;

#define EMPTY_DISTRIBUTION {0.0, 0.0, 0.0, 0.0}

inline
void
pvc_set_unscattered(pvc_distribution* dist,
		    ms_real range, ms_real rho_transmitter2)
{
  dist->power = 1.0;
  dist->power_x_width2 = range*range*rho_transmitter2;
  dist->power_x_zeta2 = rho_transmitter2;
  dist->power_x_cov = range*rho_transmitter2;
}


inline
void
pvc_add(pvc_distribution* out, pvc_distribution a, pvc_distribution b)
{
  out->power = a.power + b.power;
  out->power_x_width2 = a.power_x_width2 + b.power_x_width2;
  out->power_x_zeta2 = a.power_x_zeta2 + b.power_x_zeta2;
  out->power_x_cov2 = a.power_x_cov2 + b.power_x_cov2;
}

inline
void
pvc_subtract(pvc_distribution* out, pvc_distribution a, pvc_distribution b)
{
  out->power = a.power - b.power;
  out->power_x_width2 = a.power_x_width2 - b.power_x_width2;
  out->power_x_zeta2 = a.power_x_zeta2 - b.power_x_zeta2;
  out->power_x_cov2 = a.power_x_cov2 - b.power_x_cov2;
}

inline
void
pvc_step(pvc_distribution* dist,
	 ms_real drange, ms_real drange2, ms_real drange3,
	 ms_real ext_Theta2, ms_real trans)
{
  dist->power_x_width2 += dist->power_x_zeta2*drange2
    + 2.0 * dist->power_x_cov*drange
    + (1.0/3.0) * dist->power * ext_Theta2 * drange3;
  dist->power_x_cov += dist->power_x_zeta2*drange
    + 0.5 * dist->power * ext_Theta2 * drange;
  dist->power_x_zeta2 += ext_Theta2 * drange;
  dist->power *= trans;
  dist->power_x_width2 *= trans;
  dist->power_x_cov *= trans;
  dist->power_x_zeta2 *= trans;
}

inline
void
pvc_copy(pvc_distribution* out, pvc_distribution in)
{
  out->power = in.power;
  out->power_x_width2 = in.power_x_width2;
  out->power_x_zeta2 = in.power_x_zeta2;
  out->power_x_cov = in.power_x_cov;
}

inline
void
pvc_crop(pvc_distribution* dist, ms_real width2_max)
{
  ms_real correlation2 = dist->power_x_cov*dist->power_x_cov
    / (dist->power_x_width2*dist->power_x_zeta2);
  ms_real factor = width2_max * dist->power / dist->power_x_width2;
  ms_real factor2 = factor*factor;
  dist->power *= factor;
  dist->power_x_width2 *= factor2;
  dist->power_x_cov *= factor2;
  dist->power_x_zeta2 *= factor*(correlation2*factor + 1.0*correlation2);
}
