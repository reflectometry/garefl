#include <stdlib.h>
#include <stdio.h>
#include "refl.h"

#undef ALL_MODELS

#ifdef ALL_MODELS
  #define MODELS 5
#else
  #define MODELS 3
#endif

/*=========== CONSTRAINTS =====================*/
void constr_models(fitinfo *fit)
{
#ifdef HAVE_MAGNETIC
  int l;

  // Constraints on first data set
  fit[0].m.thetarough[1] = fit[0].m.rough[1];
  fit[0].m.thetarough[4] = fit[0].m.rough[4];
  
  fit[0].m.Prough[1]     = fit[0].m.rough[1];
  fit[0].m.Prough[2]     = fit[0].m.rough[2];
  fit[0].m.Prough[3]     = fit[0].m.rough[3];
  fit[0].m.Prough[4]     = fit[0].m.rough[4];

  fit[0].m.theta[0]      = fit[0].m.theta[1];
  fit[0].m.theta[4]      = fit[0].m.theta[3];

  // Copy constrained parameters 
  for (l=1; l < MODELS; l++) {
    fit[l].m.d[1]            = fit[0].m.d[1];
    fit[l].m.d[2]            = fit[0].m.d[2];
    fit[l].m.d[3]            = fit[0].m.d[3]; 

    fit[l].m.rho[1]        = fit[0].m.rho[1];
    fit[l].m.rho[2]        = fit[0].m.rho[2];
    fit[l].m.rho[3]        = fit[0].m.rho[3];
    fit[l].m.rho[4]        = fit[0].m.rho[4];
  
    fit[l].m.mu[1]         = fit[0].m.mu[1];
  
    fit[l].m.rough[1]      = fit[0].m.rough[1];
    fit[l].m.rough[2]      = fit[0].m.rough[2];
    fit[l].m.rough[3]      = fit[0].m.rough[3];
    fit[l].m.rough[4]      = fit[0].m.rough[4];

    fit[l].m.thetarough[1] = fit[0].m.thetarough[1];
    fit[l].m.thetarough[4] = fit[0].m.thetarough[4];
  
    fit[l].m.Prough[1]     = fit[0].m.Prough[1];
    fit[l].m.Prough[2]     = fit[0].m.Prough[2];
    fit[l].m.Prough[3]     = fit[0].m.Prough[3];
    fit[l].m.Prough[4]     = fit[0].m.Prough[4];
    
    fit[l].m.theta[0]      = fit[l].m.theta[1];
    fit[l].m.theta[4]      = fit[l].m.theta[3];
  }
  
#endif /* HAVE_MAGNETIC */
}

fitinfo* setup_models(int *models)
{
  static fitinfo fit[MODELS];
#ifdef HAVE_MAGNETIC
  int i;
  char fa[64], fb[64], fc[64], fd[64];
  fitpars *pars = &fit[0].pars;
  /* fitpars *freepars = &fit[1].pars; */
  *models = MODELS;
  char *dataset;

  for (i=0; i < MODELS; i++) fit_init(&fit[i]);
  
  /* Load the data for each model */
  
  // Second set
  dataset = "f239+245bf";   // 166G
  snprintf(fa,sizeof(fa),"%s.qa",dataset);
  snprintf(fb,sizeof(fb),"%s.qb",dataset);
  snprintf(fc,sizeof(fc),"%s.qc",dataset);
  snprintf(fd,sizeof(fd),"%s.qd",dataset);
  fit_polarized_data(&fit[0], fa, fb, fc, fd);
  
  dataset = "f257+264bf";   // 5kG
  snprintf(fa,sizeof(fa),"%s.qa",dataset);
  snprintf(fb,sizeof(fb),"%s.qb",dataset);
  snprintf(fc,sizeof(fc),"%s.qc",dataset);
  snprintf(fd,sizeof(fd),"%s.qd",dataset);
  fit_polarized_data(&fit[1], fa, fb, fc, fd);

#if 0
  dataset = "f273+279bf";   // 55G
  snprintf(fa,sizeof(fa),"%s.qa",dataset);
  snprintf(fb,sizeof(fb),"%s.qb",dataset);
  snprintf(fc,sizeof(fc),"%s.qc",dataset);
  snprintf(fd,sizeof(fd),"%s.qd",dataset);
  fit_polarized_data(&fit[2], fa, fb, fc, fd);
#endif
  dataset = "f128+134bf";   // 57G
  snprintf(fa,sizeof(fa),"%s.qa",dataset);
  snprintf(fb,sizeof(fb),"%s.qb",dataset);
  snprintf(fc,sizeof(fc),"%s.qc",dataset);
  snprintf(fd,sizeof(fd),"%s.qd",dataset);
  fit_polarized_data(&fit[2], fa, fb, fc, fd);


#ifdef ALL_MODELS
  dataset = "f167+169bf";   // 166G
  snprintf(fa,sizeof(fa),"%s.qa",dataset);
  snprintf(fb,sizeof(fb),"%s.qb",dataset);
  snprintf(fc,sizeof(fc),"%s.qc",dataset);
  snprintf(fd,sizeof(fd),"%s.qd",dataset);
  fit_polarized_data(&fit[0], fa, fb, fc, fd);

  dataset = "f109+115bf";   // 48G
  snprintf(fa,sizeof(fa),"%s.qa",dataset);
  snprintf(fb,sizeof(fb),"%s.qb",dataset);
  snprintf(fc,sizeof(fc),"%s.qc",dataset);
  snprintf(fd,sizeof(fd),"%s.qd",dataset);
  fit_polarized_data(&fit[1], fa, fb, fc, fd);

  dataset = "f159+160bf";   // 66G
  snprintf(fa,sizeof(fa),"%s.qa",dataset);
  snprintf(fb,sizeof(fb),"%s.qb",dataset);
  snprintf(fc,sizeof(fc),"%s.qc",dataset);
  snprintf(fd,sizeof(fd),"%s.qd",dataset);
  fit_polarized_data(&fit[2], fa, fb, fc, fd);

  dataset = "f128+134bf";   // 57G
  snprintf(fa,sizeof(fa),"%s.qa",dataset);
  snprintf(fb,sizeof(fb),"%s.qb",dataset);
  snprintf(fc,sizeof(fc),"%s.qc",dataset);
  snprintf(fd,sizeof(fd),"%s.qd",dataset);
  fit_polarized_data(&fit[3], fa, fb, fc, fd);

  dataset = "f49+55bf";   // 35G
  snprintf(fa,sizeof(fa),"%s.qa",dataset);
  snprintf(fb,sizeof(fb),"%s.qb",dataset);
  snprintf(fc,sizeof(fc),"%s.qc",dataset);
  snprintf(fd,sizeof(fd),"%s.qd",dataset);
  fit_polarized_data(&fit[4], fa, fb, fc, fd);
#endif

  /* Initialize instrument parameters for each model. Note
   * that if you are simultaneously fitting e.g., X-ray
   * and neutron data, then you won't want to loop here. 
   */
  for (i=0; i < MODELS; i++) {
    fit[i].beam.lambda = 4.75;
    data_resolution_fixed(&fit[i].dataA,4.75,0.021,0.,0.,0.001);
    data_resolution_fixed(&fit[i].dataB,4.75,0.021,0.,0.,0.001);
    data_resolution_fixed(&fit[i].dataC,4.75,0.021,0.,0.,0.001);
    data_resolution_fixed(&fit[i].dataD,4.75,0.021,0.,0.,0.001);
    interface_create(&fit[i].rm, "erf", erf_interface, 30);
  }

  /*============= MODEL =====================================*/

  for (i=0; i < MODELS; i++) {
  /*                         d,     rho,      mu,        rough, P,     Prough, theta, thetarough */
    //~ model_magnetic(&fit[i].m,0,     0,        0,         0,    0,        0,    0,     0);    // vacuum
    //~ model_magnetic(&fit[i].m,129.8, 7.14e-06, 1.21e-09,  9.43, 5.69e-07, 9.43, 246.8, 9.43); // crud
    //~ model_magnetic(&fit[i].m,262.4, 7.04e-06, 2.71e-09,  9.24, 1.04e-06, 9.24, 263.2, 90.0); // Fe3O4
    //~ model_magnetic(&fit[i].m,130.3, 8.25e-06, 5.73e-09,  2.48, 4.74e-06, 2.48, 272.2, 2.48); // Fe
    //~ model_magnetic(&fit[i].m,100.,  5.34e-06, 2.104e-10, 14.3, 0.,       14.3, 272.2, 14.3); // MgAl2O4

    model_magnetic(&fit[i].m,0,     0,        0,         0,    0,        0,    0,     0);    // vacuum
    model_magnetic(&fit[i].m,116.881, 7.136e-06, 2.66e-09,  9.82, 5.69e-07, 9.82, 246.8, 9.82); // crud
    model_magnetic(&fit[i].m,268.926, 6.769e-06, 2.71e-09,  43.2, 1.04e-06, 43.2, 263.2, 90.0); // Fe3O4
    model_magnetic(&fit[i].m,134.321, 8.248e-06, 5.73e-09,  11.8, 4.74e-06, 11.8, 272.2, 11.8); // Fe
    model_magnetic(&fit[i].m,100.,     5.34e-06, 2.104e-10, 3.51, 0.,       3.51, 272.2, 3.51); // MgAl2O4
  
    //~ model_magnetic(&fit[i].m,0,    0,        0,         0,    0,        0,    0,     0); // vacuum
    //~ model_magnetic(&fit[i].m,104,  7.19e-06, 1.21e-09,  5.7,  2.30e-07, 5.7,  270, 5.7); // crud
    //~ model_magnetic(&fit[i].m,280,  6.96e-06, 2.71e-09,  60.0, 8.88e-07, 60.0, 270, 60.0); // Fe3O4
    //~ model_magnetic(&fit[i].m,138,  8.27e-06, 5.73e-09,  5.4,  4.45e-06, 5.4,  270, 5.4); // Fe
    //~ model_magnetic(&fit[i].m,100., 5.10e-06, 2.104e-10, 2.4,  0.,       2.4,  270, 2.4); // MgAl2O4
  }
  
  
  /*=============== FIT PARAMETERS ===============================*/

  /* Specify which parameters are your fit parameters.
   * All parameters are placed into the fit[0] structure even if
   * they are from another model. By convention, the shared parameters
   * are all modified in model 0 and copied to the other models.
   */
  
#ifndef ALL_MODELS
  pars_add(pars,"d1",&(fit[0].m.d[1]),90.,150.);
  pars_add(pars,"d2",&(fit[0].m.d[2]),220.,300.);
  pars_add(pars,"d3",&(fit[0].m.d[3]),100.,150.);
  
  pars_add(pars,"rho2",&(fit[0].m.rho[2]),6.7e-06,7.3e-06);
  //~ pars_add(pars,"rho3",&(fit[0].m.rho[3]),7.7e-06,8.3e-06);
  pars_add(pars,"rho3",&(fit[0].m.rho[3]),7.7e-06,8.6e-06);
#endif

#ifndef ALL_MODELS
  pars_add(pars,"rho1",&(fit[0].m.rho[1]),1e-07,8.2e-06);
//  pars_add(pars,"qc1-166G2",&(fit[0].m.rho[1]),1e-07,8.2e-06);
//  pars_add(pars,"qc1-5kG",  &(fit[1].m.rho[1]),1e-07,8.2e-06);
//  pars_add(pars,"qc1-55G",  &(fit[2].m.rho[1]),1e-07,8.2e-06);
  
  pars_add(pars,"rho4",&(fit[0].m.rho[4]),5.0e-06,5.7e-06); // Should be known 5.37e-06
//  pars_add(pars,"qc4-166G2",&(fit[0].m.rho[4]),5.0e-06,5.7e-06); // Should be known 5.37e-06
//  pars_add(pars,"qc4-5kG",  &(fit[1].m.rho[4]),5.0e-06,5.7e-06); 
//  pars_add(pars,"qc4-55G",  &(fit[2].m.rho[4]),5.0e-06,5.7e-06); 
#endif

#ifdef ALL_MODELS
  pars_add(pars,"qc1-166G1",&(fit[0].m.rho[1]),1e-07,8.2e-06);
  pars_add(pars,"qc1-48G",  &(fit[1].m.rho[1]),1e-07,8.2e-06);
  pars_add(pars,"qc1-66G",  &(fit[2].m.rho[1]),1e-07,8.2e-06);
  pars_add(pars,"qc1-57G",  &(fit[3].m.rho[1]),1e-07,8.2e-06);
  pars_add(pars,"qc1-35G",  &(fit[4].m.rho[1]),1e-07,8.2e-06);
  
  pars_add(pars,"qc4-166G1",&(fit[0].m.rho[4]),5.0e-06,5.7e-06); 
  pars_add(pars,"qc4-48G",  &(fit[1].m.rho[4]),5.0e-06,5.7e-06); 
  pars_add(pars,"qc4-66G",  &(fit[2].m.rho[4]),5.0e-06,5.7e-06); 
  pars_add(pars,"qc4-57G",  &(fit[3].m.rho[4]),5.0e-06,5.7e-06); 
  pars_add(pars,"qc4-35G",  &(fit[4].m.rho[4]),5.0e-06,5.7e-06); 
#endif

  pars_add(pars,"mu1",&(fit[0].m.mu[1]),1e-10,3e-09);
  
  //~ pars_add(pars,"ro1",&(fit[0].m.rough[1]),1,15); 
  pars_add(pars,"ro1",&(fit[0].m.rough[1]),1,17); 
  pars_add(pars,"ro2",&(fit[0].m.rough[2]),1,100);
  pars_add(pars,"ro3",&(fit[0].m.rough[3]),1,15);
  pars_add(pars,"ro4",&(fit[0].m.rough[4]),1,15);  
  
#ifndef ALL_MODELS
  pars_add(pars,"th_ro2-166G2",&(fit[0].m.thetarough[2]),1,100);
  pars_add(pars,"th_ro3-166G2",&(fit[0].m.thetarough[3]),1,15);
  pars_add(pars,"th_ro2-5kG",&(fit[1].m.thetarough[2]),1,100);
  pars_add(pars,"th_ro3-5kG",&(fit[1].m.thetarough[3]),1,15);
  pars_add(pars,"th_ro2-55G",&(fit[2].m.thetarough[2]),1,100);
  pars_add(pars,"th_ro3-55G",&(fit[2].m.thetarough[3]),1,15);
#endif

#ifdef ALL_MODELS
  pars_add(pars,"th_ro2-166G1",&(fit[0].m.thetarough[2]),1,100);
  pars_add(pars,"th_ro3-166G1",&(fit[0].m.thetarough[3]),1,15);
  pars_add(pars,"th_ro2-48G",  &(fit[1].m.thetarough[2]),1,100);
  pars_add(pars,"th_ro3-48G",  &(fit[1].m.thetarough[3]),1,15);
  pars_add(pars,"th_ro2-66G",  &(fit[2].m.thetarough[2]),1,100);
  pars_add(pars,"th_ro3-66G",  &(fit[2].m.thetarough[3]),1,15);
  pars_add(pars,"th_ro2-57G",  &(fit[3].m.thetarough[2]),1,100);
  pars_add(pars,"th_ro3-57G",  &(fit[3].m.thetarough[3]),1,15);
  pars_add(pars,"th_ro2-35G",  &(fit[4].m.thetarough[2]),1,100);
  pars_add(pars,"th_ro3-35G",  &(fit[4].m.thetarough[3]),1,15);
#endif
  
#ifndef ALL_MODELS
  pars_add(pars,"qm1-166G2",&(fit[0].m.P[1]),0,5.2e-06);
  pars_add(pars,"qm2-166G2",&(fit[0].m.P[2]),0,1.5e-06);
  pars_add(pars,"qm3-166G2",&(fit[0].m.P[3]),0,6.0e-06);
  pars_add(pars,"qm1-5kG",&(fit[1].m.P[1]),0,5.2e-06);
  pars_add(pars,"qm2-5kG",&(fit[1].m.P[2]),0,1.5e-06);
  pars_add(pars,"qm3-5kG",&(fit[1].m.P[3]),0,6.0e-06);
  pars_add(pars,"qm1-55G",&(fit[2].m.P[1]),0,5.2e-06);
  pars_add(pars,"qm2-55G",&(fit[2].m.P[2]),0,1.5e-06);
  pars_add(pars,"qm3-55G",&(fit[2].m.P[3]),0,6.0e-06);
#endif

#ifdef ALL_MODELS
  pars_add(pars,"qm1-166G1",&(fit[0].m.P[1]),0,5.2e-06);
  pars_add(pars,"qm2-166G1",&(fit[0].m.P[2]),0,1.5e-06);
  pars_add(pars,"qm3-166G1",&(fit[0].m.P[3]),0,5.2e-06);
  pars_add(pars,"qm1-48G",&(fit[1].m.P[1]),0,5.2e-06);
  pars_add(pars,"qm2-48G",&(fit[1].m.P[2]),0,1.5e-06);
  pars_add(pars,"qm3-48G",&(fit[1].m.P[3]),0,5.2e-06);
  pars_add(pars,"qm1-66G",&(fit[2].m.P[1]),0,5.2e-06);
  pars_add(pars,"qm2-66G",&(fit[2].m.P[2]),0,1.5e-06);
  pars_add(pars,"qm3-66G",&(fit[2].m.P[3]),0,5.2e-06);
  pars_add(pars,"qm1-57G",&(fit[3].m.P[1]),0,5.2e-06);
  pars_add(pars,"qm2-57G",&(fit[3].m.P[2]),0,1.5e-06);
  pars_add(pars,"qm3-57G",&(fit[3].m.P[3]),0,5.2e-06);
  pars_add(pars,"qm1-35G",&(fit[4].m.P[1]),0,5.2e-06);
  pars_add(pars,"qm2-35G",&(fit[4].m.P[2]),0,1.5e-06);
  pars_add(pars,"qm3-35G",&(fit[4].m.P[3]),0,5.2e-06);
#endif
  
#ifndef ALL_MODELS
  pars_add(pars,"th1-166G2",&(fit[0].m.theta[1]),0,360);
  pars_add(pars,"th2-166G2",&(fit[0].m.theta[2]),0,360);
  pars_add(pars,"th3-166G2",&(fit[0].m.theta[3]),0,360);
  pars_add(pars,"th1-5kG",&(fit[1].m.theta[1]),0,360);
  pars_add(pars,"th2-5kG",&(fit[1].m.theta[2]),0,360);
  pars_add(pars,"th3-5kG",&(fit[1].m.theta[3]),0,360);
  pars_add(pars,"th1-55G",&(fit[2].m.theta[1]),0,360);
  pars_add(pars,"th2-55G",&(fit[2].m.theta[2]),0,360);
  pars_add(pars,"th3-55G",&(fit[2].m.theta[3]),0,360);
#endif
#ifdef ALL_MODELS
  pars_add(pars,"th1-166G1",&(fit[0].m.theta[1]),0,360);
  pars_add(pars,"th2-166G1",&(fit[0].m.theta[2]),0,360);
  pars_add(pars,"th3-166G1",&(fit[0].m.theta[3]),0,360);
  pars_add(pars,"th1-48G",&(fit[1].m.theta[1]),0,360);
  pars_add(pars,"th2-48G",&(fit[1].m.theta[2]),0,360);
  pars_add(pars,"th3-48G",&(fit[1].m.theta[3]),0,360);
  pars_add(pars,"th1-66G",&(fit[2].m.theta[1]),0,360);
  pars_add(pars,"th2-66G",&(fit[2].m.theta[2]),0,360);
  pars_add(pars,"th3-66G",&(fit[2].m.theta[3]),0,360);
  pars_add(pars,"th1-57G",&(fit[3].m.theta[1]),0,360);
  pars_add(pars,"th2-57G",&(fit[3].m.theta[2]),0,360);
  pars_add(pars,"th3-57G",&(fit[3].m.theta[3]),0,360);
  pars_add(pars,"th1-35G",&(fit[4].m.theta[1]),0,360);
  pars_add(pars,"th2-35G",&(fit[4].m.theta[2]),0,360);
  pars_add(pars,"th3-35G",&(fit[4].m.theta[3]),0,360);
#endif
  
#endif  /* HAVE_MAGNETIC */

  constraints = constr_models;

  return fit;
}
