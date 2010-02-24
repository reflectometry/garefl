/* This program is public domain. */

#ifndef _MODEL_H
#define _MODEL_H

#include "reflconfig.h"
#include "interface.h"
#include "profile.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ==== model handling ==== */

/* model_create(model*, int n)
 *   make space for n layers
 * model_extend(model*, int n)
 *   make space for n more layers
 * model_slice(model*, int l, Real d, Real rho, Real mu, Real rough)
 *   add a layer to the model
 */

#define MODEL_REPEATS_EXCEEDED -2
#define MODEL_REPEATS_BAD -3
#define MODEL_REPEATS_OVERLAP -4

#define MODEL_MAX_REPEATS 5

#ifdef HAVE_MAGNETIC
#define MODEL_FIELDS 8
#else
#define MODEL_FIELDS 4
#endif

typedef struct model_struct {
  int capacity;
  int n;
  Real *rho;
  Real *mu;
  Real *d;
  Real *rough;
#ifdef HAVE_MAGNETIC
  Real *P;
  Real *Prough;
  Real *theta;
  Real *thetarough;
#endif
  interface *rm;
  int num_repeats;
  int repeats[3*MODEL_MAX_REPEATS];
  int is_magnetic;
} model;

void model_init(model *m);
int model_extend(model *m, int n);
void model_destroy(model *m);
void model_repeat(model *m, int R, int *start, int *end, int *count);
int model_repeat_insert(model *m, int start, int end, int count);
void model_repeat_delete(model *m, int R);
void model_profile(model *m, profile *p);
void model_magnetic(model *m, Real d,
		    Real rho, Real mu, Real rough,
		    Real P, Real Prough, Real theta, Real thetarough);
void model_layer(model *m, Real d,
		 Real rho, Real mu, Real rough);
const char *model_error(int code);

void model_print(const model *m, const char *filename);

#ifdef __cplusplus
}
#endif

#endif /* _MODEL_H */

/* $Id$ */
