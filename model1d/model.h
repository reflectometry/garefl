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
 * model_slice(model*, int l, double d, double rho, double mu, double rough)
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
  double *rho;
  double *mu;
  double *d;
  double *rough;
#ifdef HAVE_MAGNETIC
  double *P;
  double *Prough;
  double *theta;
  double *thetarough;
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
void model_magnetic(model *m, double d,
		    double rho, double mu, double rough,
		    double P, double Prough, double theta, double thetarough);
void model_layer(model *m, double d, 
		 double rho, double mu, double rough);
const char *model_error(int code);

void model_print(const model *m, const char *filename);

#ifdef __cplusplus
}
#endif

#endif /* _MODEL_H */

/* $Id$ */
