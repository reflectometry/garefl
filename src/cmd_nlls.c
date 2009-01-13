/* This program is public domain */

#ifndef INCLUDING_CMD_NLLS
#error cmd_nlls.c cannot be separately compiled!
#endif

/* CVS information */
static char nlls_cvsid[] = "$Id$";

#include "nlls.h"


static int 
count_data(void)
{
  int k, nQ;

  nQ = 0;
  for (k=0; k < MODELS; k++) {
    if (fit[k].datatype == FIT_POLARIZED) 
      nQ += fit[k].dataA.n + fit[k].dataB.n + fit[k].dataC.n + fit[k].dataD.n;
    else 
      nQ += fit[k].dataA.n;
  }
  return nQ;
}

static void 
copy_all(int n, const double y[], const double dy[], double target[])
{
  int k;
  if (n==0) return;
  if (dy[0] != 0.) for (k=0; k < n; k++) target[k] = y[k]/dy[k];
  else memcpy(target, y, n*sizeof(*y));
}

static void
copy_matching(int ny, const double Qy[], const double y[],
	      int ndy, const double Qdy[], const double dy[],
	      double target[])
{
  int iy, idy, match;

  if (ny == 0 || ndy == 0) return;

  if (dy[0] != 0.) {
    match = idy = 0;
    for (iy=0; iy < ny; iy++) {
      while (idy < ndy && Qdy[idy] < Qy[iy]) idy++;
      if (Qy[iy] == Qdy[idy]) {
	target[match++] = y[iy]/dy[idy];
      }
    }
  } else {
    match = idy = 0;
    for (iy=0; iy < ny; iy++) {
      while (idy < ndy && Qdy[idy] < Qy[iy]) idy++;
      if (Qy[iy] == Qdy[idy]) {
	target[match++] = y[iy];
      }
    }
  }
}

static void 
copy_normalized_data(double *target)
{
  int k, nQ;

  nQ = 0;
  for (k=0; k < MODELS; k++) {
    if (fit[k].datatype == FIT_POLARIZED) {
      copy_all(fit[k].dataA.n,fit[k].dataA.R,fit[k].dataA.dR,target+nQ);
      nQ += fit[k].dataA.n;
      copy_all(fit[k].dataB.n,fit[k].dataB.R,fit[k].dataB.dR,target+nQ);
      nQ += fit[k].dataB.n;
      copy_all(fit[k].dataC.n,fit[k].dataC.R,fit[k].dataC.dR,target+nQ);
      nQ += fit[k].dataC.n;
      copy_all(fit[k].dataD.n,fit[k].dataD.R,fit[k].dataD.dR,target+nQ);
      nQ += fit[k].dataD.n;
    } else {
      copy_all(fit[k].dataA.n,fit[k].dataA.R,fit[k].dataA.dR,target+nQ);
      nQ += fit[k].dataA.n;
    }
  }
}

static void 
copy_normalized_theory(double *target)
{
  int k, nQ;

  nQ = 0;
  for (k=0; k < MODELS; k++) {
    if (fit[k].datatype == FIT_POLARIZED) {
      copy_matching(fit[k].nQ,fit[k].fitQ,fit[k].fitA,
		    fit[k].dataA.n,fit[k].dataA.Q,fit[k].dataA.dR,
		    target+nQ);
      nQ += fit[k].dataA.n;
      copy_matching(fit[k].nQ,fit[k].fitQ,fit[k].fitB,
		    fit[k].dataA.n,fit[k].dataB.Q,fit[k].dataB.dR,
		    target+nQ);
      nQ += fit[k].dataB.n;
      copy_matching(fit[k].nQ,fit[k].fitQ,fit[k].fitC,
		    fit[k].dataC.n,fit[k].dataC.Q,fit[k].dataC.dR,
		    target+nQ);
      nQ += fit[k].dataC.n;
      copy_matching(fit[k].nQ,fit[k].fitQ,fit[k].fitD,
		    fit[k].dataD.n,fit[k].dataD.Q,fit[k].dataD.dR,
		    target+nQ);
      nQ += fit[k].dataD.n;
    } else {
      copy_all(fit[k].dataA.n,fit[k].fitA,fit[k].dataA.dR,target+nQ);
      nQ += fit[k].dataA.n;
    }
  }
}

static void 
step_nlls(const double p[], double fx[], int n, int k, void *user_data)
{
  step_fn01(n,p,user_data);

  /* Copy normalized values into results */
  copy_normalized_theory(fx);
}

static void print_covar(int ndim, double covar[])
{
  int i,j;
  FILE *file;

  /* Print covariance matrix */
  file = fopen("covar.dat","w");
  fprintf(file,"=== Covariance matrix ===\n");
  for (i=0; i < ndim; i++) {
    for (j=0; j < ndim; j++) fprintf(file,"%12g ", covar[i*ndim+j]);
    fprintf(file,"\n");
  }
  fclose(file);
}


void cmd_nlls()
{
  double *work, *y, *bounds, *covar, *p;
  int nQ = count_data();
  int ndim = fit[0].pars.n;
  int i;

  /* Allocate storage: y, covar, bounds, p */
  work = (double *)malloc(sizeof(double)*(nQ + ndim*ndim + 3*ndim));
  assert(work != NULL);

  y = work;
  bounds = y + nQ;
  covar = bounds + 2*ndim;
  p = covar + ndim*ndim;

  write_pop(&set); /* In case fit crashes */

  /* Record what we are doing */
  if (parFD != NULL) {
    fprintf(parFD,"# %15d   Starting Levenberg-Marquardt\n", GetGen(&set));
    fflush(parFD); 
  }

  /* Suppress logging during LM */
  log_improvement = 0;

  /* Copy normalized data values */
  copy_normalized_data(y);

  /* Set bounds */
  for (i=0; i < ndim; i++) {
    bounds[i] = 0.;
    bounds[i+ndim] = 1.;
  }

  /* Get best into p */
  pars_set(&fit[0].pars, bestpars);
  pars_get01(&fit[0].pars, p);

  /* Call nlls */
  tic();
#if 1
  box_nlls(step_nlls, ndim, nQ, (double *)fit, y, bounds, p, covar);
#else
  nlls(step_nlls, ndim, nQ, (double *)fit, y, p, covar);
#endif
  printf("Done LM\n");fflush(stdout);
  toc();
  print_covar(ndim,covar);

  /* Restore logging, recording LM results */
  log_improvement = 1;
  pars_set(&fit[0].pars,bestpars);
  update_models(fit);
  record_improvement();

  /* Inject new p into the GA */
  setChromosome(&set, 0, p);

  /* Done */
  free(work);
}
