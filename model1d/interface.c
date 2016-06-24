/* This program is public domain. */

/* fzw */
/** \file
 * Handle the roughness between layers.
 */

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "interface.h"

#ifdef NEED_ERF
double erf(double);
#endif

void interface_print(interface *rm)
{
  int i;

  printf("Interface pointer %p\n",rm);
  printf("  name=%s\n",rm->name);
  printf("  n=%d\n",rm->n);
  printf("  fn=%p (tanh=%p, erf=%p)\n",rm->fn,tanh_interface,erf_interface);
  printf("  depth=%g\n",rm->depth);
  printf("  # value(%p) step(%p)\n", rm->value, rm->step);
  for (i=0; i < rm->n; i++)
    printf("  %d %g %g\n",i,rm->value[i],rm->step[i]);
}

Real interface_average(Real above, Real below, Real weight)
{
  return .5 * (above + below + (below-above) * weight);
}

Real interface_overlap(Real weight_above, Real weight_below,
				Real v_above, Real v, Real v_below)
{
  return .5 * (v_above + v
	       + (v-v_above)*weight_above
	       + v + v_below
	       + (v_below-v)*weight_below) - v;
}

/* Take derivative of z in place */
/* For first and last point use 2-point slope */
/* For other points average the slope left and right */

/* The above explanation, while correct, is not very informative.  What
 * we have coming in is the z_i value for particular f(z_i)=r_i,
 * with r_i uniformly spaced.  We want to create slices in z with these
 * z_i values near the center of the slice rather than at the edge.
 * Assuming f(z) is linear between z_i, a sensible solution is to
 * set the interfaces between the slices half way between the z_i
 * values, with the last slice extending an extra half step beyond
 * the last value.  We then convert that to a series of slice widths. */
void interface_steps(int n,Real z[])
{
   register Real ohalfstep, ztemp;
   register int j;

   ohalfstep = .5 * (z[1] - z[0]);
   for (j = 0; j < n-1; j++) {
      ztemp = z[j];
      z[j] = ohalfstep + .5 * (z[j+1] - z[j]);
      ohalfstep = .5 * (z[j+1] - ztemp);
   }
   /* Last point */
   z[j] = 2 * ohalfstep;
}

/* Comment is incorrect */
/* Constant CE that ensures that d(erf CE*Z/ZF)/dZ = .5 when Z=.5*ZF,
   where ZF is fwhm */

/* Constant CE ensures that .5*sqrt(Pi)*d(erf(y))/dy|(y=CE*Z/ZF) = .5
   when Z=.5*ZF where ZF is fwhm
   Note: .5*sqrt(Pi)*d(erf(y))/dy = exp(-y**2) */
#define CE 1.665

#if 0
/* erfinv was translated from Octave by Paul Kienzle
 * It's not the fastest or most accurate algorithm, but
 * it is good enough our purposes.
 * Copyright (C) 1995, 1996 Kurt Hornik
 * GNU General Public License 2 or later.
 */
#define TOLERANCE 1e-15
#define MAXITERATIONS 20 /* 9 should be enough */
Real erfinv(Real x)
{
  Real zero=0., s, z_old, z_new;
  int iterations;
  if (fabs(x) > 1.) return zero/zero;
  if (fabs(x) == 1.) return x/zero;
  s = sqrt(M_PI)/2;
  z_old = 1.0;
  z_new = sqrt(-log (1-fabs(x))) * copysign(1.0,x);
  iterations = 0;
  while (fabs(erf(z_new)-x) > TOLERANCE*fabs(x)) {
    z_old = z_new;
    z_new = z_old - (erf(z_old) - x) * exp(z_old*z_old) *s;
    if (++iterations > MAXITERATIONS) break;
  }
  return z_new;
}
#endif

void erf_interface(int n, Real z[], Real r[])
{
  if (n==0) {
    /* Calculate interface at a particular value */
    *r = erf(CE * *z);
  } else {
    /* Split [-1,1] into n intervals sampled at the midpoint of each. */
    Real step = 2. / (n+1);
    Real x = -1.+ step;
    int i;

    for (i=0; i<n; i++) {
#if 1
      z[i] = log((1.+x)/(1.-x)) / CE;
      r[i] = erf(z[i]*CE);
#else
      z[i] = erfinv(x)/CE;
      r[i] = x;
#endif
      x += step;
    }
    interface_steps(n,z);
  }
}

/* Comment is incorrect */
/* Ensures that d(tanh CT * Z/ ZF) / DZ = .5
   when Z = .5 * ZF, where ZF is fwhm */
/* Ensures that d(tanh(y))/dy|(y=CE*Z/ZF) = 0.334 when Z=.5*ZF,
   where ZF is fwhm */
#define CT 2.292

void tanh_interface(int n, Real z[], Real r[])
{
  if (n==0) {
    /* Calculate interface at a particular value */
    *r = tanh(CT * *z);
  } else {
    /* Split [-1,1] into n intervals sampled at the midpoint of each. */
    Real step = 2. / (n+1);
    Real x = -1.+ step;
    int i;

    for (i=0; i<n; i++) {
      z[i] = log((1.+x)/(1.-x)) / CT;
      /* FIXME isn't atanh = 0.5*(log(1.+x) - log(1.-x)) ? */
      r[i] = x;
      x += step;
    }
    interface_steps(n,z);
  }
}

/* Neville's interpolation for a quadratic. */
Real quadinterp(Real x, Real x0, Real x1, Real x2,
		  Real y0, Real y1, Real y2)
{
  Real d0 = x-x0;
  Real d1 = x-x1;
  Real d2 = x-x2;
  return (d2*(y0*d1-y1*d0)/(x0-x1) - d0*(y1*d2 - y2*d1)/(x1-x2)) / (x0-x2);
}

Real interface_value(const interface *rm, Real z)
{
  if (rm->fn != NULL) {
    /* Have an interface function */
    Real retval;
    (*(rm->fn))(0,&z,&retval);
    return retval;
  } else if (1) {
    Real retval;
    erf_interface(0,&z,&retval);
    return retval;
  } else {
    /* Need to interpolate */
    int n = rm->n;
    Real z0, z1, z2, v0, v1, v2;
    int i;
    assert(n > 1);
    /* z goes from -w/2 to w/2 in steps sampled at
     * the midpoints.  First midpoint is at:
     *     z_1 = -w/2 + w_{1}/2
     * The next midpoint is at:
     *     z_2 = z_1 + w_1
     * Curiously, if you push through the algebra you will see that
     * subsequent midpoints are at:
     *     z_{i+1}= 2 w_i + z_{i-1}
     * This makes sense considering the following profile:
     *
     *                 |yy^...
     *         |xxxx^yy|  .
     * ...^xxxx|    .  |  .
     *    .    |    .  |  .
     *    .         .     .
     *   z_{i-1}   z_i   z_{i+1}
     *
     * Here you can see that the number of xxx plus the number of yyy
     * between z_{i-1} and z_{i+1} is exactly twice xxx plus yyy, which
     * is just the width of step i.
     *
     * Finally, I'm going to put some pretend midpoints at z_1-10*w_1
     * and z_n+10*w_n with values of +/- 1 respectively.
     */
    z1 = 0.5*(rm->step[0]-rm->depth);
    z2 = z1 + rm->step[0];
    for (i=1; z > z2 && i < n; i++) {
      z0 = z1;
      z1 = z2;
      z2 = z0 + 2.*rm->step[i];
    }
    if (i==1) {
      z0 = z1 - 10.*rm->step[0];
      v0=-1.; v1=rm->value[0]; v2=rm->value[1];
    } else if (i>=n) {
      z2 = z1 + 10.*rm->step[n-1];
      v0=rm->value[n-2]; v1=rm->value[n-1]; v2=1.;
    } else {
      v0=rm->value[i-2]; v1=rm->value[i-1]; v2=rm->value[i];
    }


    if (z <= z0) return -1.0; /* way before first value */
    else if (z >= z2) return 1.0; /* way after last value */
    else return quadinterp(z,z0,z1,z2,v0,v1,v2);
  }
}

void interface_init(interface *rm)
{
  rm->step = NULL;
  rm->n = 0;
}

void interface_destroy(interface *rm)
{
  if (rm->step != NULL) free(rm->step);
  rm->step = rm->value = NULL;
  rm->n = 0;
  rm->depth = 0.0;
}

int interface_set(interface *rm, const char* name,
		  int n, const Real z[], const Real v[])
{
  rm->name = name;
  rm->n = n;
  rm->depth = 0.;
  rm->fn = NULL;
  if (n == 0) {
    /* If no roughness, then no steps needed */
    rm->step = NULL;
    return 1;
  }

  /* Step and value share one allocation */
  rm->step = malloc(sizeof(Real)*n*2);
  rm->value = rm->step+n;
  if (rm->step == NULL) {
    /* Couldn't allocate: pretend there is no roughness */
    rm->n = 0;
    return 0;
  } else {
    int i;
    /* Allocated: compute roughness steps and total depth */
    /* FIXME also want depth at 0 to support non-symmetric steps
     * and to support even numbered roughness when transitioning from
     * non-overlapping to overlapping interfaces. */
    for (i = 0; i < n; i++) {
      rm->step[i] = z[i];
      rm->value[i] = v[i];
    }
    interface_steps(n,rm->step);
    for (i = 0; i < n; i++) rm->depth += rm->step[i];
  }
  return 1;

}

int interface_create(interface *rm, const char* name,
		      interface_fn fn, int n)
{
  assert(n>=0);
  rm->name = name;
  rm->n = n;
  rm->depth = 0.;
  rm->fn = fn;
  if (n == 0) {
    /* If no roughness, then no steps needed */
    rm->step = NULL;
    return 1;
  }

  /* Step and value share one allocation */
  rm->step = malloc(sizeof(Real)*n*2);
  rm->value = rm->step+n;
  if (rm->step == NULL) {
    /* Couldn't allocate: pretend there is no roughness */
    rm->n = 0;
    return 0;
  } else {
    int i;
    /* Allocated: compute roughness steps and total depth */
    /* FIXME also want depth at 0 to support non-symmetric steps
     * and to support even numbered roughness when transitioning from
     * non-overlapping to overlapping interfaces. */
    (*fn)(n, rm->step, rm->value);
    for (i = 0; i < n; i++) rm->depth += rm->step[i];
  }
  return 1;
}

/* $Id$ */
