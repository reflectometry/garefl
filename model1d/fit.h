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
  Real lambda;         /* Wavelength */
  Real intensity;      /* Beam intensity (usually 1) */
  Real background;     /* Background signal */
  Real backabsorption; /* Absorption through substrate */
  Real Aguide;         /* Angle relative to guide */
  Real alignment;      /* Angle correction for incident beam */
  /* Note: if alignment is incorrect then the values we measure will have a
   * small non-zero Qx component.  Assuming the back slits are not too tight
   * this should not affect the total reflection measured.
   * */
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
  Real *incoherent_weights;

   /* Theory curve */
  int nQ, capacity;  /* Theory function storage */
  Real *fitQ;
  Real *fitA, *fitB, *fitC, *fitD;

  /* Reflectometry data */
  fitdata dataA, dataB, dataC, dataD;
  int datatype;
  beaminfo beam;     /* Instrument parameters */
  fitpars pars;      /* Fitting parameters */

  /* Fitness value */
  Real chisq_est;  /* Current chisq */
  Real weight;     /* Scale factor for chisq (default 1) */

	/* Working storage */
  int worksize;      /* Working storage */
  Real *work;

} fitinfo;
void fit_init(fitinfo *fit);
void fit_destroy(fitinfo *fit);
void fit_sumsq(const fitinfo *fit, int *n, Real *sumsq);
void fit_wsumsq(const fitinfo *fit, int *n, Real *sumsq);
Real fit_chisq(const fitinfo *fit);
Real fit_wchisq(const fitinfo *fit);
void fit_update(fitinfo *fit, int approx);
void fit_partial(fitinfo *fit, int approx, Real portion, Real best,
		 int weighted, int *totaldf, Real *totalsumsq);
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
typedef void fit_output(fitinfo *fit);
extern fit_output *output_model; /* constraints callback */
fitinfo*
setup_models(int *models); /* simultaneous fitting */
void
tied_parameters(fitinfo *fit); /* default constraints for simultaneous fitting */
#ifdef __cplusplus
}
#endif


#endif /* _FIT_H */

/* $Id$ */
