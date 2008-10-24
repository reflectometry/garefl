/* This program is public domain. */

#ifndef _FIT_H
#define _FIT_H

#include "model.h"
#include "data.h"
#include "profile.h"
#include "interface.h"
#include "pars.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct beaminfo_struct {
  double lambda;         /* Wavelength */
  double intensity;      /* Beam intensity (usually 1) */
  double background;     /* Background signal */
  double backabsorption; /* Absorption through substrate */
  double Aguide;         /* Angle relative to guide */
} beaminfo;

#define FIT_POLARIZED 4
#define FIT_MAGNITUDE 0
#define FIT_REAL 1
#define FIT_IMAGINARY 2

typedef struct fitinfo_struct {
	/* Reflectometry model */
  model m;           /* Model parameters */
  profile p;         /* Model profile */
  interface rm;      /* Interface shape for model */

  /* Support for incoherent sums across models */
  int number_incoherent;
  model **incoherent_models;
  double *incoherent_weights;

   /* Theory curve */
  int nQ, capacity;  /* Theory function storage */
  double *fitQ;
  double *fitA, *fitB, *fitC, *fitD;

  /* Reflectometry data */
  fitdata dataA, dataB, dataC, dataD;
  int datatype;
  beaminfo beam;     /* Instrument parameters */
  fitpars pars;      /* Fitting parameters */

  /* Fitness value */
  double chisq_est;  /* Current chisq */
  double weight;     /* Scale factor for chisq (default 1) */

	/* Working storage */  
  int worksize;      /* Working storage */
  double *work;
  
} fitinfo;
void fit_init(fitinfo *fit);
void fit_destroy(fitinfo *fit);
void fit_sumsq(const fitinfo *fit, int *n, double *sumsq);
void fit_wsumsq(const fitinfo *fit, int *n, double *sumsq);
double fit_chisq(const fitinfo *fit);
double fit_wchisq(const fitinfo *fit);
void fit_update(fitinfo *fit, int approx);
void fit_partial(fitinfo *fit, int approx, double portion, double best, 
		 int weighted, int *totaldf, double *totalsumsq);
void fit_data(fitinfo *fit, const char *file);
void fit_real_data(fitinfo *fit, const char *file);
void fit_imaginary_data(fitinfo *fit, const char *file);
void fit_polarized_data(fitinfo *fit, const char *fA, const char *fB,
		       const char *fC, const char *fD);

void fit_data_print(const fitinfo *fit);
void fit_save_staj(const fitinfo *fit, const char *file);


/* Temporary hacks to support models separate from main. */
typedef void fit_constraints(fitinfo *fit);
extern fit_constraints *constraints; /* constraints callback */
fitinfo*
setup_models(int *models); /* simultaneous fitting */
void 
tied_parameters(fitinfo *fit); /* default constraints for simultaneous fitting */
#ifdef __cplusplus
}
#endif


#endif /* _FIT_H */

/* $Id$ */
