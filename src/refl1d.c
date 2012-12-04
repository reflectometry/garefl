/* This program is public domain. */

/* Implementation and wrapping of example C fit functions */

#include <stdio.h>
#include "refl.h"


/* Application control parameters */
// FIXME: shouldn't be using globals here
//int approximate_roughness = 0;
//int weighted = 1;
//fitinfo *fit;
fit_output *output_model = NULL;
fit_constraints *constraints = NULL;
int MODELS;

/* Note: duplicated from ga_simul.c because it relies on global MODELS */
/* Copy parameters between models. */
void
tied_parameters(fitinfo fit[])
{
  int i,k;

  /* Rescue the free parameters from the model. */
  for (i=0; i < fit[1].pars.n; i++) {
    fit[1].pars.value[i] = pars_peek(&fit[1].pars,i);
  }

  /* Go through all layers copying parameters from model 0
   * to the other models. This is more work than we strictly
   * need to do (we only really need to update the parameters
   * which are varying in the fit), but the extra cost is
   * trivial compared to the cost of calculating the profile
   * and the reflectivity, so let the code be simpler.
   */
  for (k=0; k < fit[0].m.n; k++) {
    for (i=1; i < MODELS; i++) {
      fit[i].m.d[k] = fit[0].m.d[k];
      fit[i].m.rho[k] = fit[0].m.rho[k];
      fit[i].m.mu[k] = fit[0].m.mu[k];
      fit[i].m.rough[k] = fit[0].m.rough[k];
#ifdef HAVE_MAGNETIC
      fit[i].m.theta[k] = fit[0].m.theta[k];
      fit[i].m.thetarough[k] = fit[0].m.thetarough[k];
      fit[i].m.P[k] = fit[0].m.P[k];
      fit[i].m.Prough[k] = fit[0].m.Prough[k];
#endif
    }
  }

  /* Restore the free parameters to the model. */
  for (i=0; i < fit[1].pars.n; i++) {
    // Don't use pars_poke since it does bounds checking.
    *(fit[1].pars.address[i]) = fit[1].pars.value[i];
  }
}

Real
ex_update_models(fitinfo fit[], int num_models,
              int weighted, int approximate_roughness,
              int forced)
{
  int n = 0;
  Real sumsq = 0.;
  int i;
  // printf("calling update\n");
  MODELS = num_models;
  fit[0].penalty = 0.;
  if (*constraints) (*constraints)(fit);
  sumsq = fit[0].penalty;
  if (!forced && sumsq >= FIT_REJECT_PENALTY) return sumsq;
  for (i=0; i < num_models; i++) {
    int n_i = 0;
    Real sumsq_i = 0.;
    fit_update(&fit[i], approximate_roughness);
    if (weighted) fit_wsumsq(&fit[i],&n_i,&sumsq_i);
    else fit_sumsq(&fit[i],&n_i,&sumsq_i);
    fit[i].chisq_est = sumsq_i/n_i;
    n += n_i;
    sumsq += sumsq_i;
  }
  // printf("sumsq=%10g, n=%4d, pars=%d\n",sumsq,n,fit[0].pars.n);
  return n < fit[0].pars.n ? sumsq : sumsq / (n - fit[0].pars.n) ;
}


double
ex_get_penalty(fitinfo fit[])
{
  // printf("returning %.2f from penalty\n",fit[0].penalty);
  return fit[0].penalty;
}

int
ex_npars(fitinfo fit[])
{
  return fit[0].pars.n;
}

void
ex_par_values(fitinfo fit[], Real p[])
{
  pars_get(&fit[0].pars, p);
}

const char*
ex_par_name(fitinfo fit[], int i)
{
  return pars_name(&fit[0].pars, i);
}

void
ex_par_bounds(fitinfo fit[], Real lo[], Real hi[])
{
  int i;
  for (i = 0; i < fit[0].pars.n; i++) {
    lo[i] = pars_min(&fit[0].pars, i);
    hi[i] = pars_max(&fit[0].pars, i);
  }
}

void
ex_set_pars(fitinfo fit[], const Real p[])
{
  pars_set(&fit[0].pars, p);
}

int
ex_ndata(fitinfo fit[], int k, int xs)
{
  if (xs == 0)      return fit[k].dataA.n;
  else if (xs == 1) return fit[k].dataB.n;
  else if (xs == 2) return fit[k].dataC.n;
  else if (xs == 3) return fit[k].dataD.n;
  return 0;
}


const char *
ex_get_data(fitinfo fit[], int k, int xs, Real out[][4])
{
  fitdata *data;
  int i;

  if (xs == 0) data = &(fit[k].dataA);
  else if (xs == 1) data = &(fit[k].dataB);
  else if (xs == 2) data = &(fit[k].dataC);
  else if (xs == 3) data = &(fit[k].dataD);
  else return ""; // ERROR
  for (i=0; i<data->n; i++) {
    out[i][0] = data->Q[i];
    out[i][1] = data->dQ[i];
    out[i][2] = data->R[i];
    out[i][3] = data->dR[i];
  }
  return data->file;
}


void
ex_get_beam(fitinfo fit[], int k, Real values[])
{
  values[0] = fit[k].beam.lambda;
  values[1] = fit[k].beam.intensity;
  values[2] = fit[k].beam.background;
  values[3] = fit[k].beam.backabsorption;
  values[4] = fit[k].beam.Aguide;
  values[5] = fit[k].beam.alignment;
}

void
ex_set_beam(fitinfo fit[], int k, Real values[])
{
  fit[k].beam.lambda = values[0];
  fit[k].beam.intensity = values[1];
  fit[k].beam.background = values[2];
  fit[k].beam.backabsorption = values[3];
  fit[k].beam.Aguide = values[4];
  fit[k].beam.alignment = values[5];
}

int
ex_nprofile(fitinfo fit[], int k)
{
  return fit[k].p.n;
}

void
ex_get_profile(fitinfo fit[], int k,
            Real w[],
            Real rho[], Real irho[],
            Real rhoM[], Real thetaM[])
{
  int i;
  for (i = 0; i < fit[k].p.n; i++) {
    w[i] = fit[k].p.d[i];
    rho[i] = fit[k].p.rho[i];
    irho[i] = fit[k].p.mu[i]/2/fit[k].beam.lambda;
    if (fit[k].m.is_magnetic == 4) {
      rhoM[i] = fit[k].p.P[i];
      thetaM[i] = fit[k].p.theta[i];
    }
  }
}

int
ex_ncalc(fitinfo fit[], int k)
{
  return fit[k].nQ;
}

void
ex_get_reflectivity(fitinfo fit[], int k, int xs, Real Q[], Real R[])
{
  int i;
  for (i=0; i < fit[k].nQ; i++) Q[i] = fit[k].fitQ[i];
  if (xs == 0) {
    for (i=0; i < fit[k].nQ; i++) R[i] = fit[k].fitA[i];
  } else if (xs==1) {
    for (i=0; i < fit[k].nQ; i++) R[i] = fit[k].fitB[i];
  } else if (xs==2) {
    for (i=0; i < fit[k].nQ; i++) R[i] = fit[k].fitC[i];
  } else if (xs==3) {
    for (i=0; i < fit[k].nQ; i++) R[i] = fit[k].fitD[i];
  }
}

/* $Id: _cfitfn_example.cc 11 2006-03-30 08:44:57Z pkienzle $ */
