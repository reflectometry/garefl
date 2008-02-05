/* This program is public domain. */

/** \file
 *  Demonstration program for the Levenberg-Marquardt fitting routines.
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "nlls.h"

#include "randmtzig.ic"

#define BKGD 0 /* True if gaussian plus constant background */
#define BOUNDS 1 /* True if use box constraints */
#include "gaussfit.ic"


#define NDATA 20
int main(int argc, char *argv[])
{
  FILE *data;
  const int k = NDATA;
  const int n = NPAR;
  double x[2*NDATA]; /* Room for x and dy */
  double *dy = &(x[NDATA]);
  double y[NDATA], fx[NDATA];
  double p0[NPAR], p1[NPAR], p[NPAR];
  double bounds[2*NPAR];
  double covar[NPAR*NPAR];
  int i;

  if (argc > 1) init_by_int (atoi(argv[1]));
  else init_by_entropy();

  gendata(k,x,y,dy,p0);
  initp(k,x,y,p1,bounds);
  for (i=0; i < n; i++) p[i] = p1[i];
#if BOUNDS
  box_nlls(fn, n, k, x, y, bounds, p, covar);
#else
  nlls(fn, n, k, x, y, p, covar);
#endif


  printf("# ====== Pars =========================================\n");
  printf("p0: ");
  for (i=0; i < n; i++) printf("%12g ",p0[i]);
  printf("\np1: ");
  for (i=0; i < n; i++) printf("%12g ",p1[i]);
  printf("\np : ");
  for (i=0; i < n; i++) printf("%12g ",p[i]);
  printf("\nlo: ");
  for (i=0; i < n; i++) printf("%12g ",bounds[i]);
  printf("\nhi: ");
  for (i=0; i < n; i++) printf("%12g ",bounds[i+n]);
  printf("\n");
    
  printf("# ====== Covar =========================================\n");
  for (i=0; i < n; i++) {
    int j;
    for (j=0; j < n; j++) printf("%12g ",covar[i*n+j]);
    printf("\n");
  }

  printf("# ====== gnuplot command ==============================\n");
  printf("plot 'data' u 1:2:3 w errorbar, 'data' u 1:4 w l\n");

  data = fopen("data","wt");
  fprintf(data,"# x \t y \t dy \t f(x)\n");
  fn(p,fx,n,k,x);
  for (i=0; i < k; i++) {
    fprintf(data,"%g \t%g \t%g \t%g\n",x[i],y[i]*dy[i],dy[i],fx[i]*dy[i]);
  }
  fclose(data);

  exit(0);
  return 0;
}

/* $Id: test_nlls.c 35 2007-05-25 15:12:45Z ziwen $ */
