#include <stdio.h>
#include <assert.h>
#include "refl.h"
#include "reflcalc.h"
#define MODELS 2
/* initialising non-standard fitting variables */
Real vol_fract_spacer, vol_fract_alkyl, rho_spacer, rho_alkyl, global_rough;

/*=========== CONSTRAINTS =====================*/
void constr_models(fitinfo *fit)
{
  int i,k;

  /* Rescue the free parameters from the model. */
  for (i=0; i < fit[1].pars.n; i++)
    fit[1].pars.value[i] = *(fit[1].pars.address[i]);

  /* Go through all layers copying parameters from model 0 to other models */
  tied_parameters(fit);

 /* copy the global roughness to all interfaces*/
  for (i=0; i< MODELS; i++) {
    for (k=1;k<7; k++) fit[i].m.rough[k]=global_rough;
  }

  /* Restore the free parameters to the model. */
  for (i=0; i < fit[1].pars.n; i++){
    *(fit[1].pars.address[i]) = fit[1].pars.value[i];
  }

  /* allow the cr_au roughness to differ from the global rough
   * - this was saved as a free parameter for model[0], now copy to
   * model[1]
   */
  fit[1].m.rough[3]=fit[0].m.rough[3];

  /*calculate the densities for layers given in vol fractions*/
  for (i=0; i< MODELS; i++) {
    fit[i].m.rho[4]= vol_fract_spacer*rho_spacer+(1-vol_fract_spacer)*fit[i].m.rho[fit[i].m.n-1];
    fit[i].m.rho[5]= vol_fract_alkyl*rho_alkyl+(1-vol_fract_alkyl)*fit[i].m.rho[fit[i].m.n-1];
  }

  //if (fit[0].m.d[1] > 10.) {
  //  printf("fit[0].m.d[1] = %g->%g\n",fit[0].m.d[1],10.);
  //  fit[0].m.d[1] = 10.;
  //}
}

void save(fitinfo *fit)
{
  int k;

  printf("Volfrac parameters\n");
  printf("  vol_fract_spacer: %g\n",vol_fract_spacer);
  printf("  vol_fract_alkyl: %g\n",vol_fract_alkyl);
  printf("  rho_spacer: %g\n", rho_spacer);
  printf("  rho_alkyl: %g\n", rho_alkyl);
  printf("  global_rough: %g\n", global_rough);
  for (k=0; k < MODELS; k++) {
    printf("== Model %d ==\n", k);
    model_print(&fit[k].m,NULL);
  }
  
}

/*============ INITIAL SETUP ============================*/
fitinfo* setup_models(int *models)
{
  static fitinfo fit[MODELS];
  int i;
  fitpars *pars = &fit[0].pars;
  fitpars *freepars = &fit[1].pars;
  *models = MODELS;

  for (i=0; i < MODELS; i++) fit_init(&fit[i]);

  /* Load the data for each model */
  fit_data(&fit[0],"wc02.yor");
  fit_data(&fit[1],"wc03.yor");

  /* Constant resolution example */
  /* data_constant_resolution(&fit[0].dataA, 0.3); */

  /* Initialize instrument parameters for each model.*/
  for (i=0; i < MODELS; i++) {
    const Real L = 5.0042,dLoL=0.020,d=1890.0;
    Real Qlo,Tlo, dTlo,dToT,s1,s2;
    Qlo=0.0154,Tlo=0.35;
    s1=0.21, s2=s1;
    dTlo=resolution_dT(s1,s2,d);
    dToT=resolution_dToT(s1,s2,d,Tlo);
    data_resolution_fv(&fit[i].dataA,L,dLoL,Qlo,dTlo,dToT);
    fit[i].beam.lambda = L;
    fit[i].beam.background = 1e-10;
    interface_create(&fit[i].rm, "erf", erf_interface, 30);
  }

  /*============= MODEL =====================================*/

  /* Add layers: d (thickness, �), rho (Nb, �-2), mu (absorption, �?), rough (interlayer roughness, �) */
  for (i=0; i < MODELS; i++) {
    model_layer(&fit[i].m, 0.000, 2.07e-6, 0.0e-8, 10.00); /* 0 substrate */
    model_layer(&fit[i].m, 12.000, 3.60e-6, 0.0e-8, 10.00); /* 1 oxide */
    model_layer(&fit[i].m, 5.800, 3.8e-6, 0.0e-8, 10.00); /* 2 chromium */
    model_layer(&fit[i].m, 90.70, 4.50e-6, 0.0e-8, 10.00); /* 3 gold */
    model_layer(&fit[i].m, 18.00, 1.00e-6, 0.0e-8, 10.00); /* 4 spacer */
    model_layer(&fit[i].m, 28.00, -0.2e-6, 0.0e-8, 10.00); /* 5 alkyl tails */
    model_layer(&fit[i].m, 100.0, 6.35e-6, 0.0e-8, 10.00); /* 6 solvent */
  }

  rho_spacer=0.50e-6;
  rho_alkyl=-0.40e-6;
  vol_fract_spacer=0.9;
  vol_fract_alkyl=0.99;
  global_rough=10.0;

  /*correct solvent layers for different models*/
  /* fit[3].m.d[3] = ... */
  fit[1].m.rho[6]=3.4e-6;


  /*=============== FIT PARAMETERS ===============================*/

  /* Specify which parameters are your fit parameters. Parameters are fitted
   * to be the same in all datasets by default
   */
  pars_add(pars, "d_oxide", &(fit[0].m.d[1]), 8., 30.);
  pars_add(pars, "d_chromium", &(fit[0].m.d[2]), 5., 30.);
  pars_add(pars, "rho_chromium", &(fit[0].m.rho[2]), 3.0e-6, 4.5e-6);
  pars_add(pars, "rough_cr_au", &(fit[0].m.rough[3]), 5.,30.);
  pars_add(pars, "d_gold", &(fit[0].m.d[3]), 40., 150.);
  pars_add(pars, "rho_gold", &(fit[0].m.rho[3]), 4.0e-6, 4.5e-6);
  pars_add(pars, "d_spacer", &(fit[0].m.d[4]), 10., 30.);
  pars_add(pars, "rho_spacer", &(rho_spacer), 0e-6, 0.5e-6);
  pars_add(pars, "vol_fract_spacer", &(vol_fract_spacer), 0., 1.);
  pars_add(pars, "d_alkyl", &(fit[0].m.d[5]), 10., 17.);
  pars_add(pars, "rho_alkyl", &(rho_alkyl), -0.5e-6, 0.0e-6);
  pars_add(pars, "vol_fract_alkyl", &(vol_fract_alkyl), 0., 1.);
  pars_add(pars, "rho_solv_0", &(fit[0].m.rho[fit[0].m.n-1]), 6.0e-6, 6.35e-6);
  pars_add(pars, "rho_solv_1", &(fit[1].m.rho[fit[1].m.n-1]), 3.0e-6, 4.5e-6);
  pars_add(pars, "global_rough", &(global_rough), 8.0, 15.0);

  /* Build a list of 'free parameters' in fit[1].pars. These are
   * parameters for which the values are aloowed to differ from those
   * in model 0.  By default all values in all models are the same unless
   * specified here. The range data is not useful, so set it to [0,1].
   */

  pars_add(freepars, "rho_solv_1", &(fit[1].m.rho[fit[1].m.n-1]), 0., 1.);
  pars_add(freepars, "rough_cr_au", &(fit[0].m.rough[3]), 0., 1.);

  constraints = constr_models;
  output_model = save;
  return fit;
}
