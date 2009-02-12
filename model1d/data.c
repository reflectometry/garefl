/* This program is public domain. */

/* fzw */
/** \file
 * Handle data read and write.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include "data.h"
#include "reflcalc.h"

#define BAD_MEMORY -1
#define BAD_FILE -2
#define BAD_NUMBER_OF_COLUMNS -3
#define BAD_COLUMN_CONSISTENCY -4
#define BAD_RESOLUTION_VALUE -5
#define BAD_DATA_LENGTH -6
static const char* errmsg[] = {
  /*   0 */ "no error",
  /*  -1 */ "memory allocation error",
  /*  -2 */ "could not open file",
  /*  -3 */ "columns should be Q R, Q R dR, or Q dQ R dR",
  /*  -4 */ "inconsistent number of columns",
  /*  -5 */ "resolution cannot be 0",
  /*  -6 */ "no data",
} ;

const char *data_error(int code)
{
  code = -code;
  if (code < 0 || code >= sizeof(errmsg)/sizeof(errmsg[0])) code = 0;
  return errmsg[code];
}

void data_print(const fitdata *data)
{
  int i;
  for (i=0; i < data->n; i++)
    printf("%d: %g %g %g %g\n",
	   i, data->Q[i], data->dQ[i], data->R[i], data->dR[i]);
}


#define LN10 2.30258509299404568402  /* log_e 10 */
#define INF 1.e308

static void log2lin(int n, double *y, double *dy)
{
  int i;
  for (i=0; i < n; i++) {
    y[i] = exp(LN10 * y[i]);
    dy[i] = y[i] * dy[i] / LN10;
  }
}

/* If the data is substantially negative (more than 70%), then assume
   that it is log data.  If it is log data which isn't negative, then
   it is not reflectivity.  If it is linear data which is mostly negative,
   then it isn't worth fitting.
*/
int data_log2lin_if_negative(fitdata *data)
{
  int i,neg;

  neg = 0;
  for (i=0; i < data->n; i++) if (data->R[i] < 0) neg++;
  if (neg > (data->n*7)/10) {
    log2lin(data->n, data->R, data->dR);
    return 1;
  }
  return 0;
}

void data_log2lin(fitdata *data)
{
  log2lin(data->n, data->R, data->dR);
}

void data_init(fitdata *data)
{
  data->capacity = -1;
  data->n = 0;
  data->file = NULL;
}

void data_destroy(fitdata *data)
{
  data->n = 0;
  if (data->capacity > 0) {
    free(data->Q);
    data->capacity = -1;
  }
  if (data->file != NULL) {
    free(data->file);
    data->file = NULL;
  }
}

int data_create(fitdata *data, int n)
{
  if (data->capacity >= n || n==0) return 1;
  data_destroy(data);
  data->Q = malloc(4*sizeof(double)*n);
  if (data->Q != NULL) {
    data->dQ = data->Q + n;
    data->R = data->dQ + n;
    data->dR = data->R + n;
    data->capacity = n;
    return 1;
  } else {
    return 0;
  }
}

void sort_data(fitdata *data)
{
  int n;
  int data_needs_sorting = 1;

  if (data->n == 0) return;

  /* Check if the data is already sorted */
  for (n=0; n < data->n-1; n++) if (data->Q[n]>data->Q[n+1]) break;
  if (n == data->n-1) return;
  /* FIXME sort the data rather than aborting */
  data_print(data);
  assert(data_needs_sorting == 0);
}


/* Load data into data structure.
 * Return the number of data columns:
 *   2: Q, R
 *   3: Q, R, dR
 *   4: Q, dQ, R, dR
 *   dQ and dR default to 0.
 */
int data_load(fitdata *data, const char *file)
{
  FILE *f;
  int n, line, columns;
  char buf[1024];

  if (file == NULL || strlen(file)==0) return 0;

  f = fopen(file, "r");
  if (f == NULL) return BAD_FILE;

  /* Count the number of non-comment lines in the file */
  n = 0;
  while (fgets(buf,sizeof(buf),f) != NULL)
    if (buf[0] != '#') n++;

  /* Make sure there is enough room in the data block for the data */
  if (!data_create(data, n)) {
    fclose(f);
    return BAD_MEMORY;
  }

  line = columns = n = 0;
  rewind(f);
  while (fgets(buf,sizeof(buf),f) != NULL) {
    int skip = 0;
    while (isspace(buf[skip])) skip++; /* Let us skip over blank lines */
    line++;
    if (buf[skip]!='\0' && buf[skip] != '#') {
      double c1,c2,c3,c4,c5;
      int c = sscanf(buf, "%lf %lf %lf %lf %lf", &c1, &c2, &c3, &c4, &c5);
      if (c!=columns) {
	if (columns) {
	  data_destroy(data);
	  fclose(f);
	  return BAD_COLUMN_CONSISTENCY;
	} else if (c < 2 || c > 4) {
	  data_destroy(data);
	  fclose(f);
	  return BAD_NUMBER_OF_COLUMNS;
	} else {
	  columns = c;
	}
      }
      switch (c) {
      case 2:
	data->Q[n]  = c1;
	data->dQ[n] = 0.0;
	data->R[n]  = c2;
	data->dR[n] = 0.0;
	n++;
	break;
      case 3:
	if (c3 == 0.) {
	  data_destroy(data);
	  fclose(f);
	  return BAD_RESOLUTION_VALUE;
	}
	data->Q[n] = c1;
	data->dQ[n] = 0.;
	data->R[n] = c2;
	data->dR[n] = c3;
	n++;
	break;
      case 4:
	if (c2 == 0. || c4 == 0.) {
	  data_destroy(data);
	  fclose(f);
	  return BAD_RESOLUTION_VALUE;
	}
	data->Q[n] = c1;
	data->dQ[n] = c2;
	/* dQ is in sigma units. Multiply by sqrt (8 ln(2)) for FWHM */
	data->R[n] = c3;
	data->dR[n] = c4;
	n++;
	break;
      }
    }
  }
  data->n = n;
  fclose(f);
  if (n == 0) return BAD_DATA_LENGTH;

  data->file = malloc(strlen(file)+1);
  strcpy(data->file,file);
  sort_data(data);

  data->have_resolution = (columns == 4);
  data->have_uncertainty = (columns >= 3);

  /* for (n=0; n < data->n; n++) if (data->dR[n]<1e-10) data->dR[n]=1e-10; */
  return 0;
}

typedef enum {RES_VARYING, RES_FIXED, RES_DQOQ} resolution_style;
static void
calc_resolution(fitdata *data, double L, double dLoL,
		double Qlo, double Qhi, double T, resolution_style res)
{
  int lo, hi, n=data->n;
  const double *Q = data->Q;
  double *dQ = data->dQ;
  data->have_resolution = 1;
  if (n == 0) return;
  if (Q[0] < 0.) {
    if (Qhi == 0.) lo = 0;
    else for (lo = 0; lo < n && Q[lo] < -Qhi; lo++) ;
    for (hi = lo; hi < n && Q[hi] < -Qlo; hi++) ;
    /* printf("found range: [%d,%d-1] for [-%g:-%g]\n",lo,hi,Qhi,Qlo); */

    switch (res) {
    case RES_FIXED: resolution_fixed(L,dLoL,T,hi-lo,Q+lo,dQ+lo); break;
    case RES_VARYING: resolution_varying(L,dLoL,T,hi-lo,Q+lo,dQ+lo); break;
    case RES_DQOQ: resolution_dQoQ(dLoL,hi-lo,Q+lo,dQ+lo); break;
    }
  } else {
    hi = 0;
  }

  for (lo = hi; lo < n && Q[lo] < Qlo; lo++) ;
  if (Qhi == 0.) hi = n;
  else for (hi = lo; hi < n && Q[hi] < Qhi; hi++) ;
  /* printf("found range: [%d,%d-1] for [%g:%g]\n",lo,hi,Qlo,Qhi); */

  switch (res) {
  case RES_FIXED: resolution_fixed(L,dLoL,T,hi-lo,Q+lo,dQ+lo); break;
  case RES_VARYING: resolution_varying(L,dLoL,T,hi-lo,Q+lo,dQ+lo); break;
  case RES_DQOQ: resolution_dQoQ(dLoL,hi-lo,Q+lo,dQ+lo); break;
  }
}

void
data_resolution_fixed(fitdata *data, double L, double dLoL,
		      double Qlo, double Qhi, double dT)
{
  calc_resolution(data,L,dLoL,Qlo,Qhi,dT,RES_FIXED);
}

void
data_resolution_varying(fitdata *data, double L, double dLoL,
			double Qlo, double Qhi, double dToT)
{
  calc_resolution(data,L,dLoL,Qlo,Qhi,dToT,RES_VARYING);
}

void
data_resolution_fv(fitdata *data, double L, double dLoL,
		   double Qlo, double dTlo, double dToT)
{
  calc_resolution(data,L,dLoL, 0.,Qlo,dTlo,RES_FIXED);
  calc_resolution(data,L,dLoL,Qlo,0.,dToT,RES_VARYING);
}

void
data_resolution_fvf(fitdata *data, double L, double dLoL,
		    double Qlo, double Qhi,
		    double dTlo, double dToT, double dThi)
{
  calc_resolution(data,L,dLoL, 0.,Qlo,dTlo,RES_FIXED);
  calc_resolution(data,L,dLoL,Qlo,Qhi,dToT,RES_VARYING);
  calc_resolution(data,L,dLoL,Qhi, 0.,dThi,RES_FIXED);
}

void
data_resolution_dQoQ(fitdata *data, double dQoQ, double Qlo, double Qhi)
{
  calc_resolution(data,0.,dQoQ,Qlo,Qhi,0.,RES_DQOQ);
}

#if 0
/* This is old code wherein the theory line is not sampled at the same
   points as Q.  I don't know if this functionality will be needed
   later.
*/
int data_printfit_subset(const char *file, const fitdata *data,
			 int nQ, const double fitQ[], const double fitR[])
{
  FILE *f;
  int i,j;

  f = fopen(file, "w");
  fprintf(f,"# Q dQ R dR fit\n");
  if (f == NULL) return 0;
  j = 0;
  for (i=0; i < data->n; i++) {
    while (j < nQ-1 && fitQ[j] < data->Q[i]) j++;
    assert (fitQ[j] == data->Q[i]);
    fprintf(f, "%g %g %g %g %g\n",
	    data->Q[i], data->dQ[i], data->R[i], data->dR[i], fitR[j]);
  }
  fclose(f);
  return 1;
}
#endif

int data_printfit(const char *file, const fitdata *data, const double fitR[])
{
  FILE *f;
  int i;

  f = fopen(file, "w");
  fprintf(f,"# Q dQ R dR fit\n");
  if (f == NULL) return 0;
  for (i=0; i < data->n; i++) {
    fprintf(f, "%g %g %g %g %g\n",
	    data->Q[i], data->dQ[i], data->R[i], data->dR[i], fitR[i]);
  }
  fclose(f);
  return 1;
}

/* Warning: assumes fit and data are ordered from smallest to largest */
/* Warning: uses == rather than fabs(x-Tx) to align points */
void
wsumsq(const int n, const double x[], const double y[], const double dy[],
       const int Tn, const double Tx[], const double Ty[], const double Tdy[],
       int *df, double *sumsq)
{
  double v = 0.;
  int match = 0;
  int found, i, Ti;

  if (n==0 || Tn==0) return;

  if (Tn == 1) {
    /* If for some reason we are only evaluating one data point at
     * a time, we don't want to search through the entire list of
     * data points because that is an n^2 algorithm.  Instead do
     * a binary search for the first data point that matches and
     * start accumulating sumsq from there. */
    /* FIXME test case: make sure "t=wsumsq(1,Tx)" is equivalent
     * to "t=0; for (i=0; i<Tn;i++) t+=wsumsq(Tn,Tx(i))". */
    int hi=n-1, lo=0;
    while (lo < hi) {
      int mid = (lo+hi)/2;
      if (*Tx > x[mid]) lo = mid+1;
      else hi=mid;
    }
    if (x[lo] != *Tx) return;
    while (--lo >= 0 && x[lo] == *Tx) ;
    found = lo+1;
  } else {
    found = 0;
  }

  Ti = 0;
  if (Tdy != NULL) { /* Uncertainty in model and data.*/
    /* If we have uncertainty in the model but not the data then
     * something is wrong --- reverse the sense of the two in that case.
     */
    assert(dy != NULL);
    for (i=found; i < n; i++) {
      while (Ti < Tn-1 && Tx[Ti] < x[i]) Ti++;
      if (x[i] == Tx[Ti]) {
	double wdiff = (Ty[Ti]-y[i])/(Tdy[Ti]+dy[i]);
	v += wdiff*wdiff;
	match++;
      } else if (x[i] > Tx[Ti]) break;
    }
  } else if (dy != NULL) { /* Uncertainty in data. */
    for (i=found; i < n; i++) {
      while (Ti < Tn-1 && Tx[Ti] < x[i]) Ti++;
      if (x[i] == Tx[Ti]) {
	double wdiff = (Ty[Ti]-y[i])/dy[i];
	v += wdiff*wdiff;
	match++;
      } else if (x[i] > Tx[Ti]) break;
    }
  } else {  /* Uncertainty ignored. */
    for (i=found; i < n; i++) {
      while (Ti < Tn-1 && Tx[Ti] < x[i]) Ti++;
      if (x[i] == Tx[Ti]) {
	double wdiff = (Ty[Ti]-y[i]);
	v += wdiff*wdiff;
	match++;
      } else if (x[i] > Tx[Ti]) break;
    }
  }
  *df += match;
  *sumsq += v;
}

void
data_wsumsq(const fitdata *data,
	    const int fitn, const double fitQ[], const double fitR[],
	    int *n, double *sumsq)
{
  if (data->n == 0) {
    /* n && sumsq are cumulative, so no need to reset them before returning */
    return;
  }

  if (data->dR[0] == 0.) {
    /* Can't be weighted if there are no weights */
    wsumsq(data->n,data->Q,data->R,NULL,fitn,fitQ,fitR,NULL,n,sumsq);
  } else {
    wsumsq(data->n,data->Q,data->R,data->dR,fitn,fitQ,fitR,NULL,n,sumsq);
  }
}

void
data_sumsq(const fitdata *data,
	   const int fitn, const double fitQ[], const double fitR[],
	   int *n, double *sumsq)
{
  if (data->n == 0) {
    /* n && sumsq are cumulative, so no need to reset them before returning */
    return;
  }

  wsumsq(data->n,data->Q,data->R,NULL,fitn,fitQ,fitR,NULL,n,sumsq);
}

int
data_countQ(const fitdata *A, const fitdata *B,
	    const fitdata *C, const fitdata *D)
{
  int n,iA,iB,iC,iD;
  double Q, nextA, nextB, nextC, nextD;

  Q = -INF;
  n = 0;
  iA = iB = iC = iD = 0;
  nextA = A && A->n ? A->Q[0] : INF;
  nextB = B && B->n ? B->Q[0] : INF;
  nextC = C && C->n ? C->Q[0] : INF;
  nextD = D && D->n ? D->Q[0] : INF;
  while (nextA != INF || nextB != INF || nextC != INF || nextD != INF) {
    if (nextA <= nextB && nextA <= nextC && nextA <= nextD) {
      if (nextA != Q) n++, Q=nextA;
      nextA = (++iA < A->n ? A->Q[iA] : INF);
    } else if (nextB <= nextC && nextB <= nextD) {
      if (nextB != Q) n++, Q=nextB;
      nextB = (++iB < B->n ? B->Q[iB] : INF);
    } else if (nextC <= nextD) {
      if (nextC != Q) n++, Q=nextC;
      nextC = (++iC < C->n ? C->Q[iC] : INF);
    } else {
      if (nextD != Q) n++, Q=nextD;
      nextD = (++iD < D->n ? D->Q[iD] : INF);
    }
#if 0
    printf("Q=%g, A[%d]=%g, B[%d]=%g, C[%d]=%g, D[%d]=%g\n",
	   Q,iA,nextA,iB,nextB,iC,nextC,iD,nextD);
#endif
  }
  return n;
}

void
data_mergeQ(const fitdata *A, const fitdata *B,
	    const fitdata *C, const fitdata *D,
	    double merge[])
{
  int n,iA,iB,iC,iD;
  double Q, nextA, nextB, nextC, nextD;

  Q = -INF;
  n = 0;
  iA = iB = iC = iD = 0;
  nextA = A && A->n ? A->Q[0] : INF;
  nextB = B && B->n ? B->Q[0] : INF;
  nextC = C && C->n ? C->Q[0] : INF;
  nextD = D && D->n ? D->Q[0] : INF;
  while (nextA != INF || nextB != INF || nextC != INF || nextD != INF) {
    if (nextA <= nextB && nextA <= nextC && nextA <= nextD) {
      if (nextA != Q) merge[n++]=Q=nextA;
      nextA = (++iA < A->n ? A->Q[iA] : INF);
    } else if (nextB <= nextC && nextB <= nextD) {
      if (nextB != Q) merge[n++]=Q=nextB;
      nextB = (++iB < B->n ? B->Q[iB] : INF);
    } else if (nextC <= nextD) {
      if (nextC != Q) merge[n++]=Q=nextC;
      nextC = (++iC < C->n ? C->Q[iC] : INF);
    } else {
      if (nextD != Q) merge[n++]=Q=nextD;
      nextD = (++iD < D->n ? D->Q[iD] : INF);
    }
  }
}

/* $Id$ */
