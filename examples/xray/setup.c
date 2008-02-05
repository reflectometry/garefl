#include <stdio.h>
#include "refl.h"

#define MODELS 1

fitinfo* setup_models(int *models)
{
  static fitinfo fit[MODELS];
  fitpars *pars = &fit[0].pars;
  int i;

  *models = MODELS;

  for (i=0; i < MODELS; i++) fit_init(&fit[i]);

  /* Julie's data, Kevin's fit, garden path variety. */
  fit->beam.lambda = 1.54;
  fit->beam.background = 1e-10;

  fit_data(fit,"e1085009.log");
  data_log2lin(&fit[0].dataA);
  data_resolution_fixed(&fit[0].dataA,1.54,0.005,0.,0.,0.0003);
  interface_create(&fit->rm, "erf", erf_interface, 15);

  /* Add layers: depth (A), rho (Nb), mu (??), roughness (A) */
  model_layer(&fit->m, 0.0, 0.0, 0.0, 0.0); /* vacuum */
  model_layer(&fit->m, 31.8477, 8.64308823392e-05, 4.241e-05, 25.1836);
                                            /* Pt cap */
  model_layer(&fit->m, 508.784, 6.31208504302e-05, 8.241e-06, 29.9282);
                                            /* Ni80Fe20 permalloy */
  model_layer(&fit->m, 146.576, 9.38415343769e-05, 3.2246e-05, 20.2179);
                                            /* Pt55Fe45 */
  model_layer(&fit->m, 22.9417, 0.000110403595642, 4.241e-05, 20.7193);
                                            /* Seed */
  model_layer(&fit->m, 100.0, 1.50856405415e-05, 1.552e-06, 17.5265); 
                                            /* Glass */
  
  /* fit parameters */
# if 1
  pars_add(pars, "d_Pt", &(fit->m.d[1]), 0., 1000.);
  pars_add(pars, "d_NiFe", &(fit->m.d[2]), 0., 1000.);
  pars_add(pars, "d_PtFe", &(fit->m.d[3]), 0., 1000.);
  pars_add(pars, "d_Seed", &(fit->m.d[4]), 0., 1000.);
# endif

# if 1
  pars_add(pars, "Pt:Air", &(fit->m.rough[1]), 0., 50.);
  pars_add(pars, "PtFe:NiFe", &(fit->m.rough[2]), 0., 50.);
  pars_add(pars, "Seed:PtFe", &(fit->m.rough[3]), 0., 50.);
  pars_add(pars, "Glass:Seed", &(fit->m.rough[4]), 0., 50.);
# endif

# if 1
  pars_add(pars, "rho_Pt", &(fit->m.rho[1]), 0., 3e-4);
  pars_add(pars, "rho_NiFe", &(fit->m.rho[2]), 0., 3e-4);
  pars_add(pars, "rho_PtFe", &(fit->m.rho[3]), 0., 3e-4);
  pars_add(pars, "rho_Seed", &(fit->m.rho[4]), 0., 3e-4);
# endif

# if 1
  pars_add(pars, "mu_Pt", &(fit->m.mu[1]), 0., 3e-4);
  pars_add(pars, "mu_NiFe", &(fit->m.mu[2]), 0., 3e-4);
  pars_add(pars, "mu_PtFe", &(fit->m.mu[3]), 0., 3e-4);
  pars_add(pars, "mu_Seed", &(fit->m.mu[4]), 0., 3e-4);
# endif

  return fit;
}
