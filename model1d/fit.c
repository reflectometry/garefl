/* This program is public domain. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "refl.h"
#include "reflcalc.h"

static double frandom(void) { return (double)rand()/RAND_MAX; }

/** \file
 * Correct for intensity, background and for Q<0, absorption through
 * the substrate.  Back absorption must be applied before the resolution
 * calculation.  Constant intensity and background can be applied before
 * or after resolution with the same result, so for convenience all are
 * done within the same function.
 *
 * If intensity were to vary from measurement to measurement then it
 * would have to be applied after the resolution calculation.  Similarly
 * for the portion of background which varies from measurement to 
 * measurement independent of incident angle.  The portion of background
 * due to properties of the sample such as off-specular scattering and 
 * sample warp are angle dependent and should be applied before the 
 * resolution calculation.
 */
void 
beam_apply(const beaminfo *beam, int N, const double Q[], double R[])
{
  const double r=beam->backabsorption;
  const double I=beam->intensity;
  const double B=beam->background;
  int i;

  for (i=0; i < N; i++) R[i] = (Q[i]<0?r:1.) * I * R[i] + B; 
}

void 
beam_init(beaminfo *beam)
{
  beam->backabsorption = 1.;
  beam->intensity = 1.;
  beam->background = 0.;
  beam->lambda = 1.;
  beam->Aguide = -90.;
}

void 
fit_init(fitinfo *fit)
{
  beam_init(&fit->beam);
  pars_init(&fit->pars);
  profile_init(&fit->p);
  model_init(&fit->m);
  interface_init(&fit->rm);
  fit->m.rm = &fit->rm;
  fit->capacity = -1;
  fit->nQ = 0;
  data_init(&fit->dataA);
  data_init(&fit->dataB);
  data_init(&fit->dataC);
  data_init(&fit->dataD);
  fit->worksize = 0;
  fit->datatype = FIT_MAGNITUDE;
  fit->weight = 1.;

  /* Parameters to support incoherent sum of models */
  fit->number_incoherent = 0;
  fit->incoherent_models = NULL;
  fit->incoherent_weights = NULL;
}

static void save_gj2(FILE *fid, const fitinfo *fit)
{
#ifdef HAVE_MAGNETIC
  int i;

  /* Don't support repeats yet */
  assert(fit->m.num_repeats == 0);
  fprintf(fid, "%#15.6lE%#15.6lE%#15.6lE\n", fit->beam.lambda, 0.01, 0.0005);
                                           /* Made up numbers for dL and dT */
  fprintf(fid, "%#15.6E%#15.6E\n",fit->beam.intensity, fit->beam.background);
  fprintf(fid, "%13d%13d%13d\n", fit->m.n-1, fit->m.rm->n,0);
  fprintf(fid, "%#15.6E%#15.6E%13d\n", 0., 0.5, 200); /* Qa range */
  fprintf(fid, "%#15.6E%#15.6E%13d\n", 0., 0.5, 200); /* Qb range */
  fprintf(fid, "%#15.6E%#15.6E%13d\n", 0., 0.5, 200); /* Qc range */
  fprintf(fid, "%#15.6E%#15.6E%13d\n", 0., 0.5, 200); /* Qd range */
  fprintf(fid, "%s\n", "E"); /* roughness profile */
  fprintf(fid, "%s\n", " abcd"); /* cross sections */
  fprintf(fid, "%.*s\n", (int)strlen(fit->dataA.file)-1, fit->dataA.file);
  fprintf(fid, "%s\n", ""); /* outfile */

  for (i = 0; i < fit->m.n; i++) {
    fprintf(fid, "%#15.6lE%#15.6lE%#15.6lE%#15.6lE\n",
	    fit->m.rho[i]*16.*M_PI,fit->m.d[i], fit->m.rough[i], fit->m.mu[i]);
    fprintf(fid, "%#15.6lE%#15.6lE%#15.6lE\n",
	    fit->m.P[i]*16.*M_PI, fit->m.d[i], fit->m.Prough[i]);
    fprintf(fid, "%#15.6lE\n",
	    fit->m.theta[i]);
  }

  fprintf(fid, "\n\n");
#endif
}

static void save_mlayer_vacuum(FILE *fid,int n)
{
  int i;
  for (i=0; i < n; i++)
    fprintf(fid, "%#15.6lE%#15.6lE%#15.6lE%#15.6lE%#15.6lE\n",
	    0.,0.,0.,0.,0.);
}

static void save_mlayer_layers(FILE *fid, const model *m, int start, int stop)
{
  int i;
  for (i = start; i < stop; i++)
    fprintf(fid, "%#15.6lE%#15.6lE%#15.6lE%#15.6lE%#15.6lE\n",
	    m->rho[i]*16.*M_PI,0.,m->d[i],m->rough[i],m->mu[i]);
}

static void save_mlayer(FILE *fid, const fitinfo *fit)
{
  int layers, top, mid, bottom;

  /* Don't support repeats yet */
  assert(fit->m.num_repeats == 0);
  assert(fit->m.n <= 28);

  /* Separate into vacuum-top-mid-bottom */
  layers = fit->m.n;
  top=mid=bottom=1;
  layers -= 4;
  if (layers > 0) {
    if (layers <= 8) {
      mid += layers;
    } else if (layers <= 16) {
      mid = 9;
      top += layers-8;
    } else {
      mid = top = 9;
      bottom += layers-16;
    }
  }

  /* Write header */
  fprintf(fid,"%13d%13d%13d%13d%13d%13d\n",
	  top, mid, bottom, 1, 0, fit->m.rm->n);
  fprintf(fid, "%#15.6G%#15.6G%#15.6G\n",
	  fit->beam.lambda, 0.01, 0.0005); /* Made up numbers for dL and dT */
  fprintf(fid, "%#15.6G%#15.6G%#15.6G%#15.6G%15d\n",
	  fit->beam.intensity, fit->beam.background, 
	  0., 0.5, 200); /* Qmin Qmax #Q */
  fprintf(fid, "%s\n", "e");
  fprintf(fid, "%s\n", fit->dataA.file);
  fprintf(fid, "%s\n", "");

  /* Write layers */
  layers = fit->m.n;
  if (layers == 1) {
    save_mlayer_vacuum(fid,5);            /* vacuum up to substrate */
    save_mlayer_layers(fid,&(fit->m),0,1); /* substrate */
  } else if (layers == 2) {
    save_mlayer_layers(fid,&(fit->m),0,1); /* Vacuum */
    save_mlayer_layers(fid,&(fit->m),0,1); /* Top == vacuum */
    save_mlayer_vacuum(fid,1);            /* Top-mid break */
    save_mlayer_layers(fid,&(fit->m),1,2); /* Mid == substrate */
    save_mlayer_vacuum(fid,1);            /* Mid-bottom break */
    save_mlayer_layers(fid,&(fit->m),1,2); /* Bottom == substrate */
  } else if (layers == 3) {
    save_mlayer_layers(fid,&(fit->m),0,1); /* Vacuum */
    save_mlayer_layers(fid,&(fit->m),0,1); /* Top == vacuum */
    save_mlayer_vacuum(fid,1);            /* Top-mid break */
    save_mlayer_layers(fid,&(fit->m),1,2); /* Mid */
    save_mlayer_vacuum(fid,1);            /* Mid-bottom break */
    save_mlayer_layers(fid,&(fit->m),2,3); /* Bottom == substrate */
  } else {
    save_mlayer_layers(fid,&(fit->m),0,1+top);          /* Vacuum+top layers*/
    save_mlayer_vacuum(fid,1);                         /* Top-mid break */
    save_mlayer_layers(fid,&(fit->m),1+top,1+top+mid);  /* Mid layers */
    save_mlayer_vacuum(fid,1);                         /* Mid-bottom break */
    save_mlayer_layers(fid,&(fit->m),1+top+mid,layers); /* Bottom layers */
  }
    
  /* No fit parameters */
  fprintf(fid, "\n");
}

void fit_save_staj(const fitinfo *fit, const char *filename)
{
  FILE *fid;

  if (fit->m.n == 0) return; /* Quietly abort if no layers */

  if (filename == NULL) {
    fid = stdout;
  } else {
    fid = fopen(filename, "w");
    if (fid == NULL) return;
  }

  if (fit->datatype == FIT_POLARIZED) save_gj2(fid,fit);
  else save_mlayer(fid,fit);

  if (fid != stdout) fclose(fid);
}



void fit_destroy(fitinfo *fit)
{
  pars_destroy(&fit->pars);
  profile_destroy(&fit->p);
  model_destroy(&fit->m);
  data_destroy(&fit->dataA);
  data_destroy(&fit->dataB);
  data_destroy(&fit->dataC);
  data_destroy(&fit->dataD);
  interface_destroy(&fit->rm);
  if (fit->worksize>0) free(fit->work);
  if (fit->capacity>0) free(fit->fitQ);
  fit->capacity = -1;
  fit->nQ = 0;
  fit->worksize = 0;
}

void fit_wsumsq(const fitinfo *fit, int *n, double *sumsq)
{
  if (fit->datatype == FIT_POLARIZED) {
    data_wsumsq(&fit->dataA, fit->dataA.n, fit->dataA.Q, fit->fitA, n, sumsq);
    data_wsumsq(&fit->dataB, fit->dataB.n, fit->dataB.Q, fit->fitB, n, sumsq);
    data_wsumsq(&fit->dataC, fit->dataC.n, fit->dataC.Q, fit->fitC, n, sumsq);
    data_wsumsq(&fit->dataD, fit->dataD.n, fit->dataD.Q, fit->fitD, n, sumsq);
  } else {
#undef DOTEST
#ifdef DOTEST
    int i,n1=0,n2=0;
    double sum1=0.,sum2=0.;
    /* Test case to check that single point sumsq is working. */
    printf("Testing single point sumsq\n");
    for (i=0; i < fit->dataA.n; i++)
      data_wsumsq(&fit->dataA, 1, fit->dataA.Q+i, fit->fitA+i, n1, sum1);
    data_wsumsq(&fit->dataA, fit->dataA.n, fit->dataA.Q, fit->fitA, n2, sum2);
    assert(n1==n2);
    assert(sum1==sum2);
    n += n1;
    sumsq += sum1;
#else
    data_wsumsq(&fit->dataA, fit->dataA.n, fit->dataA.Q, fit->fitA, n, sumsq);
#endif
  }
  *sumsq *= fit->weight;
}

void fit_sumsq(const fitinfo *fit, int *n, double *sumsq)
{
  if (fit->datatype == FIT_POLARIZED) {
    data_sumsq(&fit->dataA, fit->dataA.n, fit->dataA.Q, fit->fitA, n, sumsq);
    data_sumsq(&fit->dataB, fit->dataB.n, fit->dataB.Q, fit->fitB, n, sumsq);
    data_sumsq(&fit->dataC, fit->dataC.n, fit->dataC.Q, fit->fitC, n, sumsq);
    data_sumsq(&fit->dataD, fit->dataD.n, fit->dataD.Q, fit->fitD, n, sumsq);
  } else {
    data_sumsq(&fit->dataA, fit->dataA.n, fit->dataA.Q, fit->fitA, n, sumsq);
  }
  *sumsq *= fit->weight;
}

double fit_chisq(const fitinfo *fit)
{
  double sumsq = 0.;
  int n = 0;
  fit_sumsq(fit,&n,&sumsq);
  return sumsq / (n - fit->pars.n);
}

double fit_wchisq(const fitinfo *fit)
{
  double sumsq = 0.;
  int n = 0;
  fit_wsumsq(fit,&n,&sumsq);
  return sumsq / (n - fit->pars.n);
}


static void getdata(fitdata *data, const char *file)
{
  int c = data_load(data,file);
  if (c < 0) {
    /* FIXME don't print errors and abort from the library */
    /* maybe record that there is an error in the fit structure so */
    /* that subsequent calls to fit are ignored even if user ignores */
    /* the return value from the data_load procedure. */
    fprintf(stderr,"%s: %s\n",file,data_error(c));
    exit(1);
  }
}

void
fit_polarized_data(fitinfo *fit, const char *fileA, const char *fileB, 
		   const char *fileC, const char *fileD)
{
  getdata(&fit->dataA,fileA);  
  getdata(&fit->dataB,fileB);  
  getdata(&fit->dataC,fileC);  
  getdata(&fit->dataD,fileD);  
  fit->nQ = -1;
  fit->datatype = FIT_POLARIZED;
}

void
fit_data(fitinfo* fit, const char *file)
{
  getdata(&fit->dataA,file);
  fit->nQ = -1;
  fit->datatype = FIT_MAGNITUDE;
}

void
fit_real_data(fitinfo* fit, const char *file)
{
  getdata(&fit->dataA,file);
  fit->nQ = -1;
  fit->datatype = FIT_REAL;
}

void
fit_imaginary_data(fitinfo* fit, const char *file)
{
  getdata(&fit->dataA,file);
  fit->nQ = -1;
  fit->datatype = FIT_IMAGINARY;
}

void fit_data_print(const fitinfo *fit)
{
  if (fit->datatype == FIT_POLARIZED) {
    printf("== A ==\n");
    data_print(&fit->dataA);
    printf("== B ==\n");
    data_print(&fit->dataB);
    printf("== C ==\n");
    data_print(&fit->dataC);
    printf("== D ==\n");
    data_print(&fit->dataD);
  } else {
    data_print(&fit->dataA);
  }
}

/* Find Q values needed to calculate the fit, including padding at the
 * ends to help calculate the resolution accurately.  For thick layers,
 * need to insert Q points in the middle as well. 
 */
static void find_target_Q(fitinfo *fit)
{
  int nQ;

  if (fit->nQ >= 0) return;
  if (fit->datatype == FIT_POLARIZED) {
    nQ = data_countQ(&fit->dataA,&fit->dataB,&fit->dataC,&fit->dataD);
    if (5*nQ >= fit->capacity) {
      if (fit->capacity > 0) free(fit->fitQ);
      fit->fitQ = malloc(5*nQ*sizeof(double));
      assert(fit->fitQ != NULL); /* FIXME don't abort */
      fit->capacity = 5*nQ;
      fit->fitA = fit->fitQ+nQ;
      fit->fitB = fit->fitQ+2*nQ;
      fit->fitC = fit->fitQ+3*nQ;
      fit->fitD = fit->fitQ+4*nQ;
    }
    fit->nQ = nQ;
    data_mergeQ(&fit->dataA,&fit->dataB,&fit->dataC,&fit->dataD,fit->fitQ);
  } else {
    /* Need countQ/mergeQ to remove duplicates with same Q, different dQ */
    nQ = data_countQ(&fit->dataA,NULL,NULL,NULL);
    if (2*nQ >= fit->capacity) {
      if (fit->nQ) free(fit->fitQ);
      fit->fitQ = malloc(2*nQ*sizeof(double)); 
      assert(fit->fitQ != NULL);
      fit->capacity = 2*nQ;
      fit->fitA = fit->fitQ+nQ;
      fit->fitB = fit->fitC = fit->fitD = NULL;
    }
    fit->nQ = nQ;
    data_mergeQ(&fit->dataA,NULL,NULL,NULL,fit->fitQ);
  }
  /*{ int i; for (i=0; i < nQ; i++) printf("Q[%d]: %g\n",i,fit->fitQ[i]); }*/
}

static void apply_beam_parameters(fitinfo *fit, double *cross_section,
				  fitdata *data)
{
  /* Account for back absorption, background and intensity.
   * Strictly speaking background and intensity should be
   * applied after resolution is calculated, but mathematically
   * it is equivalent to apply them before or after.  Back absorption,
   * which varies with Q, must be applied before resolution.
   */
  beam_apply(&fit->beam, fit->nQ, fit->fitQ, cross_section);

  /* Compute the convolution of the reflectivity */
  if (fit->work != NULL) {
    /* FIXME We are not allocating enough storage for multiple measurements
       with different resolution.  Make sure that we fail in a controlled
       way until we have a chance to fix this. */
    assert(data->n <= fit->nQ);
    resolution(fit->nQ, fit->fitQ, cross_section,
               data->n, data->Q, data->dQ, fit->work);
    memcpy(cross_section, fit->work, sizeof(double)*data->n);
  } else {
    assert(data->n == fit->nQ);
    /* FIXME --- with no resolution, need to copy/interpolate raw data
       values at data->Q from the computed values at fit->fitQ.
    */
  }
  /* FIXME sumsq is now based on data Q values rather than
   * theory Q values, so we have grounding problems. */
  assert(fit->work != NULL);
}

static void extend_work(fitinfo *fit, int worksize)
{
  if (fit->worksize < worksize) {
    if (fit->worksize) free(fit->work);
    fit->work = malloc(sizeof(double)*worksize);
    assert(fit->work != NULL);
    fit->worksize = worksize;
  }
}

static void calc_approx(fitinfo *fit)
{
  /* Don't have approximate roughness model for polarized beam yet */
  assert(fit->datatype == FIT_MAGNITUDE);

  /* Generate reflectivity amplitude from the profile */
  reflrough(fit->m.n, fit->m.d, fit->m.rough, fit->m.rho, 
	    fit->m.mu, fit->beam.lambda,
	    fit->nQ, fit->fitQ, fit->fitA);
  apply_beam_parameters(fit,fit->fitA,&fit->dataA);
}


/* For debugging, writes profile to file */
void _write_profile(const fitinfo *fit, const char file[])
{
  int i;
  FILE *f = fopen(file,"w");
  assert(f!=NULL);

  if (fit->datatype == FIT_POLARIZED) {
#ifdef HAVE_MAGNETIC
    fprintf(f,"\n# %12s %12s %12s %12s %12s %12s\n",
	    "depth","rho","mu","P","theta","cos(theta");
    for (i=0; i < fit->p.n; i++)
      fprintf(f,"%12.7f %12.7f %12.7f %12.7f %12.7f %12.7f\n",
	      fit->p.d[i],fit->p.rho[i],fit->p.mu[i],fit->p.P[i],
	      fit->p.theta[i],fit->p.expth[2*i]);
#endif
  } else {
    fprintf(f,"\n# %12s %12s %12s\n",
	    "depth","rho","mu");
    for (i=0; i < fit->p.n; i++)
      fprintf(f,"%12.7f %12.7f %12.7f\n",
	      fit->p.d[i],fit->p.rho[i],fit->p.mu[i]);
  }
  fclose(f);
}

void _write_refl(const fitinfo *fit, const char name[])
{
  int i,j;
  FILE *f = fopen(name,"w");
  assert(f!=NULL);

  if (fit->datatype == FIT_POLARIZED) {
    fprintf(f,"\n# %10s %10s %18s %18s %18s %18s\n","Q","dQ","A","B","C","D");
    j = 0;
    for (i=0; i < fit->nQ; i++) {
      double dQ;
      while (j<fit->dataA.n && fit->dataA.Q[j] < fit->fitQ[i]) j++;
      dQ = (j<fit->dataA.n && fit->dataA.Q[j] == fit->fitQ[i] ? 
	    fit->dataA.dQ[j]: 0.);
      fprintf(f, "%10.5f %10.5f %18.5g %18.15g %18.15g %18.15g\n", 
	      fit->fitQ[i],dQ,
	      fit->fitA[i],fit->fitB[i],fit->fitC[i],fit->fitD[i]);
    }
  } else {
    fprintf(f,"\n# %10s %10s %18s\n","Q","dQ","R");
    j = 0;
    for (i=0; i < fit->nQ; i++) {
      double dQ;
      while (j<fit->dataA.n && fit->dataA.Q[j] < fit->fitQ[i]) j++;
      dQ = (j<fit->dataA.n && fit->dataA.Q[j] == fit->fitQ[i] ? 
	    fit->dataA.dQ[j]: 0.);
      fprintf(f, "%10.5f %10.5f %18.5g\n", fit->fitQ[i],dQ,fit->fitA[i]);
    }
  }
  fclose(f);
}

/* Incoherent sum of multiple models for magnetic reflectometry */
static void incoherent_magnetic_theory(fitinfo *fit)
{
	double total_weight; /* Total weight of all models */
	int i, k;
	profile p; /* Incoherent model profile */
	double *A,*B,*C,*D;

	/* Make space for incoherent models. */
 	extend_work(fit,4*fit->nQ);
 	A = fit->work;
	B = fit->work + fit->nQ;
	C = fit->work + 2*fit->nQ;
	D = fit->work + 3*fit->nQ;

	/* Incoherent sum of the theory functions. */
  profile_init(&p);
	for (i=0; i < fit->number_incoherent; i++) {
	  model_profile(fit->incoherent_models[i], &p);
	  magnetic_reflectivity(p.n, p.d, p.rho, p.mu, fit->beam.lambda,
		  p.P, p.expth, fit->beam.Aguide, fit->nQ, fit->fitQ, A,B,C,D);
	  for (k=0; k < fit->nQ; k++) {
	  	fit->fitA[k] += A[k] * fit->incoherent_weights[i];
	  	fit->fitB[k] += B[k] * fit->incoherent_weights[i];
	  	fit->fitC[k] += C[k] * fit->incoherent_weights[i];
	  	fit->fitD[k] += D[k] * fit->incoherent_weights[i];
	  }
	}
  profile_destroy(&p);

	/* Incoherent sum of models requires relative model weighting.
	 * The base model is assumed to have weight 1.  The remaining
	 * models can have their weights adjusted, with the result
	 * normalized by the total weight. */
	total_weight = 1.;
	for (i=0; i < fit->number_incoherent; i++) {
		total_weight += fit->incoherent_weights[i];
	}
	for (k=0; k < fit->nQ; k++) {
		fit->fitA[k] /= total_weight;
		fit->fitB[k] /= total_weight;
		fit->fitC[k] /= total_weight;
		fit->fitD[k] /= total_weight;
	}
}

/* Incoherent sum of multiple models for magnetic reflectometry */
static void incoherent_nonmagnetic_theory(fitinfo *fit)
{
	double total_weight; /* Total weight of all models */
	int i, k;
	profile p; /* Incoherent model profile */
	double *A;

	/* Make space for incoherent models. */
 	extend_work(fit,fit->nQ);
 	A = fit->work;

	/* Incoherent sum of the theory functions. */
  profile_init(&p);
	for (i=0; i < fit->number_incoherent; i++) {
		profile_reset(&p);
	  model_profile(fit->incoherent_models[i], &p);
	  /* profile_print(&p,NULL); */
	  reflectivity(p.n, p.d, p.rho, p.mu, fit->beam.lambda,
		  fit->nQ, fit->fitQ, A);
	  for (k=0; k < fit->nQ; k++) {
	  	fit->fitA[k] += A[k] * fit->incoherent_weights[i];
	  }
	}
  profile_destroy(&p);

	/* Incoherent sum of models requires relative model weighting.
	 * The base model is assumed to have weight 1.  The remaining
	 * models can have their weights adjusted, with the result
	 * normalized by the total weight. */
	total_weight = 1.;
	for (i=0; i < fit->number_incoherent; i++) {
		total_weight += fit->incoherent_weights[i];
	}
	for (k=0; k < fit->nQ; k++) {
		fit->fitA[k] /= total_weight;
	}
}

static void calc_magnitude(fitinfo *fit)
{
  /* _write_profile(fit,"prof.out"); */

  if (fit->datatype == FIT_POLARIZED) {
#ifdef HAVE_MAGNETIC
    /* We've got polarized data: use magnetic calculations */
    magnetic_reflectivity(fit->p.n, fit->p.d, fit->p.rho, 
			  fit->p.mu,fit->beam.lambda,
			  fit->p.P, fit->p.expth, 
			  fit->beam.Aguide, fit->nQ, fit->fitQ,  
			  fit->fitA, fit->fitB, fit->fitC, fit->fitD);

	  if (fit->number_incoherent > 0) incoherent_magnetic_theory(fit);

    /* _write_refl(fit,"reflA.out"); */
    apply_beam_parameters(fit,fit->fitA,&fit->dataA);
    apply_beam_parameters(fit,fit->fitB,&fit->dataB);
    apply_beam_parameters(fit,fit->fitC,&fit->dataC);
    apply_beam_parameters(fit,fit->fitD,&fit->dataD);
    /* _write_refl(fit,"reflB.out"); */
#else  /* !HAVE_MAGNETIC */
    fprintf(stderr,"Refllib was not compiled for magnetic systems.\n");
    exit(1);
#endif /* !HAVE_MAGNETIC */
  } else {
    /* Generate reflectivity amplitude from the profile */
    reflectivity(fit->p.n, fit->p.d, fit->p.rho, 
				fit->p.mu, fit->beam.lambda,
				fit->nQ, fit->fitQ, fit->fitA);

	  if (fit->number_incoherent > 0) incoherent_nonmagnetic_theory(fit);

    /* _write_refl(fit,"reflA.out"); */
    apply_beam_parameters(fit,fit->fitA,&fit->dataA);
    /* _write_refl(fit,"reflB.out"); */
  }

}

static void calc_real(fitinfo *fit)
{
  assert(fit->datatype == FIT_REAL);

  /* Generate reflectivity amplitude from the profile */
  reflectivity_real(fit->p.n, fit->p.d, fit->p.rho, 
		    fit->p.mu, fit->beam.lambda,
		    fit->nQ, fit->fitQ, fit->fitA);	

  /* FIXME need to support repeated Q with different resolution. For
     the resolution case, this simply means making sure that enough
   */
  assert(fit->nQ == fit->dataA.n); 

  /* Compute the convolution of the reflectivity */
  /* FIXME can you simply convolve the real reflectivity? */
  if (fit->work != NULL) {
    resolution(fit->nQ, fit->fitQ, fit->fitA,
               fit->dataA.n, fit->dataA.Q, fit->dataA.dQ, fit->work);
    memcpy(fit->fitA, fit->work, sizeof(double)*fit->dataA.n);
  }

  /* FIXME sumsq is now based on data Q values rather than
   * theory Q values, so we have grounding problems. */
  assert(fit->work != NULL);
}


static void calc_imaginary(fitinfo *fit)
{
  assert(fit->datatype == FIT_IMAGINARY);

  /* Generate reflectivity amplitude from the profile */
  reflectivity_imag(fit->p.n, fit->p.d, fit->p.rho, 
		    fit->p.mu, fit->beam.lambda,
		    fit->nQ, fit->fitQ, fit->fitA);

  /* FIXME need to support repeated Q with different resolution. For
     the resolution case, this simply means making sure that enough
   */
  assert(fit->nQ == fit->dataA.n); 

  /* Compute the convolution of the reflectivity */
  /* FIXME can you simply convolve the real reflectivity? */
  if (fit->work != NULL) {
    resolution(fit->nQ, fit->fitQ, fit->fitA,
               fit->dataA.n, fit->dataA.Q, fit->dataA.dQ, fit->work);
    memcpy(fit->fitA, fit->work, sizeof(double)*fit->dataA.n);
  }

  /* FIXME sumsq is now based on data Q values rather than
   * theory Q values, so we have grounding problems. */
  assert(fit->work != NULL);
}

#if 0 /* Experimental code to evaluate only a portion of the Q values */
/* Strongly favour low Q just beyond Qc */
static double Qweight(double minusQc, double plusQc, double Q)
{
  if (Q <= minusQc) return 1./(1+minusQc-Q);
  else if (Q < plusQc) return 0.01;
  else return 1./(1+Q-plusQc);
}

/* Generate approximately the correct portion of Q values */
static int
generate_subset(const int n, const double Q[], double subQ[],
		const double mQc, const double pQc, const double portion)
{
  int i, j;
  double total = 0.;

  /* Compute probability for each Q */
  for (i=0; i < n; i++) total += Qweight(mQc,pQc,Q[i]);
  total *= portion;

  /* Determine which Q to choose */
  j=0;
  for (i=0; i < n; i++)
    if (Qweight(mQc,pQc,Q[i])/total<frandom()) subQ[j++] = Q[i];

  /* Return number of Q chosen */
  return j;
}

void fit_portion_update(fitinfo *fit, double portion)
{
  int worksize = 0;

  /* Generate the profile from the model */
  model_profile(&fit->m, &fit->p);

  /* For convolution, we need one double per Q value. */
  if (worksize < fit->totalQ) worksize += fit->totalQ;
  extend_work(fit,worksize);

  if (portion > 0.9) { /* Must save at least 10% */
    /* Cheat: don't copy subQ.  Instead need to be able to find
     * the memory allocated for it, which is allQ+n.
     */
    fit->fitQ = fit->allQ;
    fit->nQ = fit->totalQ;
  } else {
    /* Weight Q selection according to Qc */
    double rhoV = fit->m.rho[0];
    double rhoS = fit->m.rho[fit->m.n-1];
    double pQc, mQc;
    if (rhoS<=rhoV) {
      mQc = -sqrt(16.*M_PI*(rhoV-rhoS));
      pQc = 0.;
    } else {
      mQc = 0.;
      pQc = sqrt(16.*M_PI*(rhoS-rhoV));
    }
    fit->fitQ = fit->allQ + fit->totalQ;
    fit->nQ = generate_subset(fit->totalQ, fit->allQ ,fit->fitQ, 
			      mQc, pQc, portion);
  }

  if (fit->datatype == FIT_REAL) calc_real(fit);
  else if (fit->datatype == FIT_IMAGINARY) calc_imaginary(fit);
  else calc_magnitude(fit);
}
#endif

void fit_update(fitinfo *fit, int approx)
{
  int worksize = 0;

  /* Find the Q points at which we need to calculate the theory */
  find_target_Q(fit);

  /* For convolution, we need one double per Q value. */
  if (worksize < fit->nQ) worksize += fit->nQ;
  extend_work(fit,worksize);

  /* Calculate the theory */
  model_profile(&fit->m, &fit->p);
  if (fit->datatype == FIT_REAL) calc_real(fit);
  else if (fit->datatype == FIT_IMAGINARY) calc_imaginary(fit);
  else if (approx && fit->datatype != FIT_POLARIZED) calc_approx(fit);
  else calc_magnitude(fit);
}

#define NOVALUE 1e308
static void 
partial_point(fitinfo *fit, int k, double Q[],
	      double A[], double B[], double C[], double D[],
	      int weighted, int *df, double *sumsq)
{
  /* FIXME don't remove this without fixing fit_w?sumsq.  In particular,
     it is now assuming fitQ matches data.Q, which is not true for this
     code.  The underlying problem is that the theory calculation points
     not correspond to the data points, either because there are multiple
     points with the same resolution or because all four cross sections
     are not measured at every Q value or because thick layers requires
     oversampling in Q to avoid aliasing effects.
   */
  assert(1==0); 
  if (A[k] == NOVALUE) {
    /* Set next point */
    fit->fitQ = Q+k;
    fit->fitA = A+k;
    fit->fitB = B+k;
    fit->fitC = C+k;
    fit->fitD = D+k;
    
    /* Compute reflectivity */
    if (fit->datatype == FIT_REAL) calc_real(fit);
    else if (fit->datatype == FIT_IMAGINARY) calc_imaginary(fit);
    else calc_magnitude(fit);
    
    /* Accumulate sumsq */
    if (weighted) fit_wsumsq(fit,df,sumsq);
    else fit_sumsq(fit,df,sumsq);
  }
}

void
fit_partial(fitinfo *fit, int approx, double portion, double best, 
	    int weighted, int *totaldf, double *totalsumsq)
{
  int i, samples;
  double *work, *Q, *A, *B, *C, *D;
  int nQ, worksize;
  double sumsq = 0.;
  int df = 0;

  /* FIXME kill partial for now --- we need to refactor the
   * the other bits so they allow single point computations a little
   * more cleanly.
   */
  portion = 1.;
  /* If all, update the model and calculate sumsq using all points. */
  if (portion >= 1.) {
    fit_update(fit, approx);
    if (weighted) fit_wsumsq(fit,&df,&sumsq);
    else fit_sumsq(fit,&df,&sumsq);
    fit->chisq_est = sumsq/df;
    return;
  }

  /* Generate the profile from the model */
  model_profile(&fit->m, &fit->p);

  /* Find the Q points at which we need to calculate the model */
  find_target_Q(fit);

  /* Clear work vector so resolution isn't calculated */
  worksize = fit->worksize; fit->worksize = 0;
  work = fit->work; fit->work = NULL;

  /* Stash Q,R info */
  nQ = fit->nQ;
  Q = fit->fitQ;
  A = fit->fitA;
  B = fit->fitB;
  C = fit->fitC;
  D = fit->fitD;

  /* For the rest of this function, assume a single point Q vector */
  fit->nQ = 1;

  /* Clear results vector */
  for (i=0; i < nQ; i++) A[i] = NOVALUE;

  /* First pass --- calculate a few randomly selected Q values */
  /* We wat to take the minimum number of samples which will allow us to */
  /* reject the hypothesis that the current individual comes from a */
  /* chisq distribution as good or better than the one which gave rise */
  /* to the best chisq.   Any individual so bad is unlikely to propagate */
  /* its genes to the next generation. For individuals as good as or */
  /* better than the best, we want to keep sampling up to some portion. */
  /* FIXME use a proper statistical test? */
  for (i=0; i < 5; i++) {
    int k= (int)floor(frandom()*nQ);
    partial_point(fit,k,Q,A,B,C,D,weighted,&df,&sumsq);
  }
  samples = ceil(portion*nQ);
  while (i++ < samples) {
    int k= (int)floor(frandom()*nQ);
    partial_point(fit,k,Q,A,B,C,D,weighted,&df,&sumsq);
    if (sumsq/df > 10.*best) break;
  }

  /* Accumulate the sumsq for this fit. */
  /* printf("sample %d of %d with %g vs. %g\n",df,samples,sumsq/df,best); */
  *totalsumsq += sumsq;
  *totaldf += df;
  fit->chisq_est = sumsq/df;

  /* Restore work vector */
  fit->worksize = worksize;
  fit->work = work;

  /* Restore Q,R info */
  fit->nQ = nQ;
  fit->fitQ = Q;
  fit->fitA = A;
  fit->fitB = B;
  fit->fitC = C;
  fit->fitD = D;
}

/* $Id$ */
