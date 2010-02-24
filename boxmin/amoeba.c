/*
 * This program is public domain.
 */

/** \file  
  * Nelder-Mead Downhill Simplex Algorithm, use same scale for all contractions.
	*
	* Based on amoeba from Numerical Recipes Second Edition 
	* by Press, Tekolsky, Vetterling and Flannery
	*/

#undef DEBUG
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "amoeba.h"

#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

/* Support for noisy fn using vertex averaging */
#define AVG 1
#if AVG>1
#define MEASURE(s,p) (avg(s->fn,s->n,p,s->userdata))
static Real avg(optimfn fn, int ndim, Real p[], void *data)
{
  int j;
  Real v = 0.;

  for (j=0; j<AVG; j++) v += (*fn)(ndim,p,data);
  return v/AVG;
}
#else
#define MEASURE(s,p) ((*s->fn)(s->n,p,s->userdata))
#endif

static void metropolitan_init(simplex *s, Real scale)
{
  const Real *bounds = s->bounds;
  Real *po = amoeba_VERTEX(s,0);
  int ndim = s->n;
  int i;
  Real v;

  /* Duplicate P_o in each vertex of simplex */
  for (i=1; i <= ndim; i++) 
	memcpy(amoeba_VERTEX(s,i),po,ndim*sizeof(Real));
  
  /* Step out by a factor of scale */
  for (i=0; i < ndim; i++) {
    if (bounds == NULL) {
      /* Take a step along dim i */
      v = po[i]+scale;
    } else {
      /* Determine step size from scale and bounds */
      Real lo = bounds[i], hi=bounds[i+ndim];
      Real step = scale*(hi-lo)/100.;

      /* Take a step along dim i, keeping it in bounds */
      v = po[i];
      v += step;
      if (v > hi) {
	v -= 2*step;
	if (v < lo) v = lo;
      }
    }
    amoeba_VERTEX(s,i+1)[i] = v;
  }

  /* Measure objective function at each vertex */
  for (i=0; i <= ndim; i++) 
    amoeba_VALUE(s,i) = MEASURE(s, amoeba_VERTEX(s,i));
}


static void sumvertices(simplex *s)
{
  Real *psum = s->psum;
  int ndim = s->n;
  int i, j;
  Real *pi;

  for (j=0; j<ndim; j++) psum[j] = 0.;
  for (i=0; i<=ndim;i++) {
    pi = amoeba_VERTEX(s,i);
    for (j=0; j<ndim; j++) psum[j] += pi[j];
  }
}

static Real trymove(simplex *s, int ihi, Real scale)
{
  int ndim = s->n;
  const Real *bounds = s->bounds;
  Real *psum = s->psum;
  Real *ptry = s->ptry;
  Real *phi = amoeba_VERTEX(s,ihi);

  Real fac1, fac2, ytry;
  int j;


  /* We are moving the worst vertex in the simplex in the direction */
  /* of the center of the simplex (excluding the worst vertex) by a */
  /* factor of scale (which may be negative if we are moving away). */
  /* Given that psum contains the sum of the vertices, the */
  /* following code does this, though it is not obvious from */
  /* looking at it. */
  fac1 = (1.-scale)/ndim;
  fac2 = fac1-scale;
  if (bounds == NULL) {
    for (j=0; j<ndim; j++) ptry[j] = psum[j]*fac1 - phi[j]*fac2;
  } else {
    for (j=0; j<ndim; j++) {
      ptry[j] = psum[j]*fac1 - phi[j]*fac2;
      if (ptry[j] < bounds[j]) ptry[j] = bounds[j];
      else if (ptry[j] > bounds[j+ndim]) ptry[j] = bounds[j+ndim];
    }
  }

  /* If the new vertex is an improvement, use it. */
  ytry = MEASURE(s,ptry);
  if (ytry < amoeba_VALUE(s,ihi)) {
    amoeba_VALUE(s,ihi) = ytry;
    for (j=0; j<ndim; j++) {
      psum[j] += ptry[j] - phi[j]; /* update sum substituting new hi */
      phi[j] = ptry[j]; /* replace hi with new hi */
    }
  }
  return ytry;
}

/* Find worst (ihi), next worst (inhi) and best (ilo) points on simplex */
static void find_extremes(simplex *s)
{
  int i, ihi,inhi,ilo;
  int ndim = s->n;
  
  ilo=0;
  ihi = amoeba_VALUE(s,0)>amoeba_VALUE(s,1) ? 0 : 1;
  inhi = 1-ihi;
  for (i=0;i<=ndim;i++) {
    Real yi = amoeba_VALUE(s,i);
    if (yi < amoeba_VALUE(s,ilo)) ilo=i;
    if (yi > amoeba_VALUE(s,ihi)) {
      inhi=ihi;
      ihi=i;
    } else if (yi > amoeba_VALUE(s,inhi)) {
      if (i != ihi) inhi=i;
    }
  }
  s->ilo = ilo; s->ihi = ihi; s->inhi = inhi;
}

void amoeba_dumpsimplex(simplex *s)
{
  int i, j;
  Real *po = amoeba_VERTEX(s,0);
  Real sumsq = 0.;

  printf("---Simplex---\n");
  for (i=0; i <= s->n; i++) {
    Real *pi = amoeba_VERTEX(s,i);
    printf("P[%d]=%g at ( ",i,amoeba_VALUE(s,i));
    for (j=0; j<s->n; j++) {
      printf("%g ",pi[j]);
      sumsq += (po[j]-pi[j])*(po[j]-pi[j]);
    }
    printf(")\n");
  }
  printf("volume ~ %g\n",sqrt(sumsq));
}

/* returns n^2+4n+1 */
int amoeba_worksize(int n)
{
  return n*n + 4*n + 1;
}

void amoeba_init
(simplex *s,
 int n,                 /* Problem dimension */
 const Real bounds[], /* 2*n bounds as lo 1,2,3 ... hi 1 2 3 ...,or NULL */
 Real work[],         /* min. n^2+4n+1 values for simplex */
 optimfn fn,            /* Function to optimize */
 void *userdata)        /* User data structure */
{
  s->iterations = 0;
  s->fn = fn;
  s->userdata = userdata;
  s->n = n;
  s->bounds = bounds;    
  s->p = work; /* Simplex must be at the beginning */
  s->ptry = s->p + (n+1)*(n+1);
  s->psum = s->ptry + n;
}

void amoeba_reset(simplex *s, Real Po[], Real scale)
{
  int i;

  /* Put the starting value into the first simplex */
  for (i=0; i < s->n; i++) s->p[i] = Po[i];
  /* Put the starting value into the first simplex */
  metropolitan_init(s,scale);
  sumvertices(s);
}


Real amoeba_flatness(simplex *s)
{
  Real yhi = amoeba_VALUE(s,s->ihi);
  Real ylo = amoeba_VALUE(s,s->ilo);

  return (yhi == ylo) ? 0. : 2. * fabs(yhi-ylo) / (fabs(yhi)+fabs(ylo));
}


Real* amoeba_best(simplex *s)
{
  find_extremes(s);
  return amoeba_VERTEX(s,s->ilo);
}

Real* amoeba_step(simplex *s)
{
  int i, j;
  int ilo = s->ilo;
  int ihi = s->ihi;
  int inhi = s->inhi;
  int ndim = s->n;
  Real ytry;
  Real ylo = amoeba_VALUE(s,ilo);
  Real yhi = amoeba_VALUE(s,ihi);
  Real ynhi = amoeba_VALUE(s,inhi);

#ifdef DEBUG
    amoeba_dumpsimplex(s);
#endif

  /* Transform simplex */
  /* printf("alpha transform\n"); */
  ytry = trymove(s,ihi,-ALPHA);
  
  if(ytry < ylo) {
    /* printf("best...extrapolate\n"); */
    trymove(s,ihi,GAMMA);
  } else if (ytry >= ynhi) {
    /* printf("worse...try closer\n"); */
    ytry = trymove(s,ihi,BETA);
    if (ytry >= yhi) {
      /* printf("closer is not better...suck in all vertices\n"); */
      Real *plo = amoeba_VERTEX(s,ilo);
      
#if 0
      Real size = 0.;
      Real *phi = amoeba_VERTEX(s,ihi);
      for(j=0;j<ndim;j++) size += (phi[j]-plo[j])*(phi[j]-plo[j]);
      printf("size = %g\n",size);
      if (size < ftol*ftol) {
	printf("simplex size < %g\n", ftol);
	return ilo;
      }
#endif
      
      for(i=0;i<=ndim;i++) {
	if(i != ilo) {
	  Real *ptmp = amoeba_VERTEX(s,i);
	  /* original code: p[i] = (p[i] + p[ilo])*0.5; */
	  /* Norm Berk mod: p[i] = p[i]*beta + p[ilo]*(1.0-beta); */
	  for (j=0;j<ndim;j++) ptmp[j] = ptmp[j]*BETA + plo[j]*(1.-BETA);
	  amoeba_VALUE(s,i) = MEASURE(s,ptmp);
	}
      }
      sumvertices(s);
    }
  }

  return amoeba_best(s);
}

Real *amoeba(simplex *s,       /* Simplex */
	       Real ftol,      /* Minimum flatness of simplex */
	       int itmax)        /* Maximum iterations */
{
  int iterations = 0;
  Real *best;

  best = amoeba_best(s);
  for(;;) {
    if (amoeba_flatness(s)<ftol) break;
    best = amoeba_step(s);
    if (iterations++ >= itmax) break;
  }

  return best;
}

/* $Id: amoeba.c 35 2007-05-25 15:12:45Z ziwen $ */
