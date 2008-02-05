/* This program is public domain. */

/** \file
 *	Demonstration program for the amoeba fitting routines.
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "amoeba.h"

#include "randmtzig.ic"
#include "peaks.ic"

int main(int argc, char *argv[])
{
  int i,k;
  double p[NDIM*NDIM+4*NDIM+1];
  double dup[(NDIM+1)*(NDIM+1)];
  double bounds[NDIM*2];
  simplex S;
  simplex *s = &S;
  double *pbest;

  srand(time(NULL));
  init_by_entropy();

#if 0
  amoeba_init(s,NDIM,NULL,p,fobjective,NULL);
#else
  amoeba_init(s,NDIM,bounds,p,fobjective,NULL);
#endif

  /* random initialization within bounds */
  setbounds(NDIM,bounds);
  for (k=0; k <= NDIM; k++) {
    double *pk = amoeba_VERTEX(s,k);
    for (i=0; i < NDIM; i++) {
      double pos = (double)(rand()%10000)/10000.0;
      pk[i] = pos*(bounds[i+NDIM]-bounds[i]) + bounds[i];
    }
  }

#if 0
  /* Near global minimum */
  p[0] = MIN_X0 + 0.005;
  p[1] = MIN_Y0 + 0.005;
#endif

  /* Calculate the function values at the vertices */
  fncalls = 0;
  for (i=0; i <= NDIM; i++) 
    amoeba_VALUE(s,i) = fobjective(NDIM,amoeba_VERTEX(s,i),NULL);

  /* Save a copy of the initial simplex; normally not needed */
  memcpy(dup,p,(NDIM+1)*(NDIM+1)*sizeof(double));

  amoeba_dumpsimplex(s);
  pbest = amoeba_best(s);
  for (i=0; i < 1000; i++) {
    if (amoeba_flatness(s) < 1e-10) break;
    pbest = amoeba_step(s);
    if (i%100 == 99) amoeba_dumpsimplex(s);
  }

  printf("================================================\n");
  printf("\ninitial simplex\n");
  s->p=dup; amoeba_dumpsimplex(s); s->p=p; /* temporarily swap out simplex */
  printf("\nfinal simplex\n");
  amoeba_dumpsimplex(s);

  printf("\n\ncalls=%d\n",fncalls);
  printf("best=%g found at (%g",pbest[NDIM],pbest[0]);
  for (i=1; i < NDIM; i++) printf(",%g",pbest[i]);
  printf(")\n");
  targets();
  neighbourhood("bowl",NDIM,pbest,fobjective,NULL);
  exit(0);
  return 0;
}

/* $Id: test_amoeba.c 35 2007-05-25 15:12:45Z ziwen $ */
