/* This program is public domain. */

#include "reflconfig.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "profile.h"

#include <stdio.h>
void profile_print(profile *p, char *filename)
{
  int i, step;
  FILE *f;
  if (filename == NULL) {
    f = stdout;
    step = 0;
  } else {
    f = fopen(filename, "w");
    if (f == NULL) return;
    step = 1;
  }

  if (step) {
    Real z = -p->vacuum_offset;
#ifdef HAVE_MAGNETIC
    fprintf(f,"# z rho mu P theta\n");
    for (i=0; i < p->n; i++) {
      fprintf(f, "%g %g %g %g %g\n",z,p->rho[i],p->mu[i],p->P[i],p->theta[i]);
      z += p->d[i];
      fprintf(f, "%g %g %g %g %g\n",z,p->rho[i],p->mu[i],p->P[i],p->theta[i]);
    }
#else
    fprintf(f,"# z rho mu\n");
    for (i=0; i < p->n; i++) {
      fprintf(f, "%g %g %g\n",z,p->rho[i], p->mu[i]);
      z += p->d[i];
      fprintf(f, "%g %g %g\n",z,p->rho[i], p->mu[i]);
    }
#endif
  } else {
#ifdef HAVE_MAGNETIC
    fprintf(f,"# z rho mu P theta\n");
    for (i=0; i < p->n; i++)
      fprintf(f, "%d %g %g %g %g %g\n",i,p->d[i],
	      p->rho[i],p->mu[i],p->P[i],p->theta[i]);
#else
    fprintf(f,"# z rho mu\n");
    for (i=0; i < p->n; i++)
      fprintf(f, "  %d %g %g %g\n",i,p->d[i],p->rho[i],p->mu[i]);
#endif
  }
  if (f != stdout) fclose(f);

}

#define PROFILE_DEFAULT_LENGTH 1000

void profile_init(profile *p)
{
  p->capacity = p->n = 0;
  p->vacuum_offset = 0.;
}

int profile_good(profile *p) { return (p->capacity >= 0); }

void profile_slice(profile *p, Real d, Real rho, Real mu
#ifdef HAVE_MAGNETIC
		   , Real P, Real theta
#endif
		   )
{
  int n = p->n;
#ifdef VERBOSE
#ifdef HAVE_MAGNETIC
  printf("%d: %g %g %g %g %g\n", n, d, rho, mu, P, theta);
#else
  printf("%d: %g %g %g\n", n, d, rho, mu);
#endif
#endif
#ifdef HAVE_MAGNETIC
  p->P[n] = P;
  p->theta[n] = theta;
#endif
  p->rho[n] = rho;
  p->mu[n] = mu;
  p->d[n] = d;
  p->n++;
}

void profile_expth(profile *p)
{
#ifdef HAVE_MAGNETIC
  int n = p->n;
  int i;

  for (i=0; i < n; i++) {
    const Real th = p->theta[i];
    if (th != 270.) {
# if defined(HAVE_SINCOS) && !defined(USE_SINGLE)
      sincos(th*M_PI/180.,p->expth+2*i+1,p->expth+2*i);
# else
      p->expth[2*i] = cos(th*M_PI/180.);
      p->expth[2*i+1] = sin(th*M_PI/180.);
# endif
    } else {
      p->expth[2*i] = 0.;
      p->expth[2*i+1] = -1.;
    }
  }
#endif
}


/* Reset the profile, recreating it if it was bad. */
void profile_reset(profile *p)
{
  if (p->capacity < 0) p->capacity = 0;
  p->n = 0;
  p->vacuum_offset = 0.;
}

/* Clear all profile memory and reinitialize the structure. */
void profile_destroy(profile *p)
{
  if (p->capacity > 0) {
    if (p->mu != NULL) free(p->mu);
  }
  profile_init(p);
}

/* Extend the size of a profile to accomodate new slices. */
int profile_extend(profile *p, int n)
{
  int i;

  if (p->capacity < 0) return 0;
  if (p->n+n > p->capacity) {
    int old_size = p->capacity;
    int new_size = p->capacity + n;
    new_size += new_size/10; /* 10% spare */
    if (new_size < PROFILE_DEFAULT_LENGTH) new_size = PROFILE_DEFAULT_LENGTH;
    if (p->capacity <= 0) {
      p->mu = malloc(PROFILE_FIELDS*sizeof(Real)*new_size);
    } else {
      p->mu = realloc(p->mu,PROFILE_FIELDS*sizeof(Real)*new_size);
      for (i=PROFILE_FIELDS-1; i>0; i--)
	memmove(p->mu+i*new_size,p->mu+i*old_size,old_size*sizeof(Real));
    }
    p->rho = p->mu + new_size;
    p->d = p->mu + 2*new_size;
#ifdef HAVE_MAGNETIC
    p->P = p->mu + 3*new_size;
    p->theta = p->mu + 4*new_size;
    p->expth = p->mu + 5*new_size;
#endif
    p->capacity = new_size;
    if (p->mu==NULL) {
      profile_destroy(p);
      p->capacity = -1; /* remember that we destroyed the profile */
      return 0;
    }
  }
  return 1;
}

/* total profile depth, including vacuum interface */
Real profile_depth(profile *p)
{
  int i;
  Real z=0.;
  for (i=0; i < p->n; i++) z += p->d[i];
  return z;
}

/* Repeated sections: copy the individual profile slices
 * from one part of the profile to another.
 */
void profile_copy(profile *p, int start, int length)
{
  int i, end=p->n;
  if (!profile_extend(p, length)) return;
  for (i=0; i < PROFILE_FIELDS; i++) {
    Real *field = p->mu + i*p->capacity;
    memcpy(field+end, field+start, sizeof(Real)*length);
  }
  p->n += length;
}

/* $Id$ */
