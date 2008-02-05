/* This program is public domain. */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "pars.h"

/// CVS information
static char cvsid[] = "$Id$";

void pars_init(fitpars *pars)
{
  pars->capacity = 0;
  pars->n = 0;
}

void pars_reset(fitpars *pars)
{
  if (pars->capacity < 0) pars->capacity = 0;
  pars->n = 0;
}

int pars_extend(fitpars *pars, int n)
{
  int size = pars->capacity;
  int need = pars->n + n;

  if (size < 0) return 0;
  if (need > size) {
    /* Since we have to allocate anyway, add 10% spare plus some */
    need += need/10 + 20;
    /* Allocate or extend existing memory */
    if (size <= 0) {
      pars->min = malloc(3*sizeof(double)*need);
      pars->address = malloc(2*sizeof(void*)*need);
    } else {
      pars->address = realloc(pars->address, 2*sizeof(void*)*need);
      pars->min = realloc(pars->min, 3*sizeof(double)*need);
      /* Shift columns to new offsets, starting at the last one. */
      if (pars->address != NULL) {
	memmove(pars->address+need, pars->address+size, size*sizeof(void*));
      }
      if (pars->min != NULL) {
	memmove(pars->min+2*need, pars->min+2*size, size*sizeof(double));
	memmove(pars->min+need, pars->min+size, size*sizeof(double));
      }
    }
    /* Check that memory was allocated successfully */
    if (pars->min==NULL || pars->address==NULL) {
      pars_destroy(pars);
      pars->capacity = -1; /* remember that we destroyed */
      return 0;
    }
    /* Update column offsets */
    pars->capacity = need;
    pars->range = pars->min + need;
    pars->value = pars->range + need;
    pars->name = (const char **)(pars->address) + need;
  }
  return 1;
}

void pars_destroy(fitpars *pars)
{
  if (pars->capacity > 0) {
    if (pars->min != NULL) free(pars->min);
    if (pars->address != NULL) free(pars->address);
  }
  pars_init(pars);  
}

double* pars_vector(fitpars *pars)
{
  return pars->value;
}
double pars_to_01(const fitpars *pars, int i, double v)
{
  double min=pars->min[i];
  if (v < min) return 0.;
  else if (v > min+pars->range[i]) return 1.;
  else return (v-min)/pars->range[i];
}
double pars_from_01(const fitpars *pars, int i, double v)
{
  return pars->min[i]+v*pars->range[i];
}
const char* pars_name(const fitpars *pars, int i)
{
  return pars->name[i];
}
double pars_peek(const fitpars *pars, int i)
{
  return *(pars->address[i]);
}
double pars_peek01(const fitpars *pars, int i)
{
  return pars_to_01(pars,i,*(pars->address[i]));
}
void pars_poke(const fitpars *pars, int i, double v)
{
  double min = pars->min[i];
  double max = min + pars->range[i];
  pars->value[i] = *(pars->address[i]) = v < min ? min : (v > max ? max : v);
}
void pars_poke01(const fitpars *pars, int i, double v)
{
  pars->value[i] = (v < 0. ? 0. : (v > 1. ? 1. : v));
  *(pars->address[i]) = pars_from_01(pars,i,pars->value[i]);
}
double pars_min(const fitpars *pars, int i)
{
  return pars->min[i];
}
double pars_max(const fitpars *pars, int i)
{
  return pars->min[i]+pars->range[i];
}
void pars_set(const fitpars *pars, const double v[])
{
  int i;
  for (i=0; i < pars->n; i++) pars_poke(pars,i,v[i]);
}
void pars_set01(const fitpars *pars, const double v[])
{
  int i;
  for (i=0; i < pars->n; i++) pars_poke01(pars,i,v[i]);
}
void pars_get(const fitpars *pars, double v[])
{
  int i;
  for (i=0; i < pars->n; i++) v[i] = pars_peek(pars,i);
}
void pars_get01(const fitpars *pars, double v[])
{
  int i;
  for (i=0; i < pars->n; i++) v[i] = pars_peek01(pars,i);
}
int pars_count(const fitpars *pars)
{
  return pars->n;
}


void pars_print_set(fitpars *pars, double *set)
{
  int i;
  for (i=0; i<pars->n; i++) {
    char range[]="..........";
    double value, portion;
    if (pars->mode == PARS_ZERO_ONE) {
      portion = set[i];
      value = set[i]*pars->range[i]+pars->min[i];
    } else {
      portion = (set[i]-pars->min[i])/pars->range[i];
      value = set[i];
    }
    if (portion < 0.) portion = 0.; 
    else if (portion >= 1.) portion = 0.999999;

    printf("%3d ",i);
    if (pars->name[i] != NULL) printf("%25s", pars->name[i]);
    range[(int)floor(portion*strlen(range))] = '|';
    printf(" %s ", range);
    printf("%g in [%g,%g]\n",value,pars->min[i],pars->min[i]+pars->range[i]);
  }
}

void pars_print(fitpars *pars)
{
  pars_print_set(pars, pars->value);
}

int pars_select(fitpars *pars)
{
  char buffer[100];
  int p;

  while (1) {
    printf("Enter parameter number: "); fflush(stdout);
    if (fgets(buffer,sizeof(buffer),stdin) == NULL) return -1;
    if (sscanf(buffer,"%d",&p) == 1) {
      if (p >= 0 && p < pars->n) return p;
    }
    if (isalpha(buffer[0])) return -1;
    pars_print(pars);
  }
}

int pars_enter_value(fitpars *pars, int p, double *v)
{
  char buffer[100];

  if (p >= 0) {
    const char *name = pars->name[p];
    double min = pars->min[p];
    double max = min + pars->range[p];
    float val;
    printf("Enter value for %s in [%g,%g]: ",name,min,max);
    fflush(stdout);
    if (fgets(buffer,sizeof(buffer),stdin)!=NULL) {
      if (sscanf(buffer,"%f",&val)>0 && val>=min && val<=max) {
	*v = val;
      } else {
	p = -1;
	printf("Invalid value ignored.\n");
      }
    }
  }
  return p;
}


void pars_add(fitpars *pars, const char *name, double *d, 
	      double min, double max)
{
#if 0
  printf("add par %d: ",pars->n);
  if (name != NULL) printf("%s=", name, name);
  printf("%g in [%g,%g]\n",*d,min,max);
#endif

  pars_extend(pars,1);
  if (pars->n < pars->capacity) {
    pars->name[pars->n] = name;
    pars->address[pars->n] = d;
    pars->min[pars->n] = min;
    pars->range[pars->n] = max-min;
    pars->n++;
  }
}

void pars_set_zero_one(fitpars *pars)
{
  int i;

  /* Translate real numbers in the model to parameters */
  pars->mode = PARS_ZERO_ONE;
  for (i=0; i < pars->n; i++) {
    pars->value[i] = pars_to_01(pars, i, *(pars->address[i]));
  }
}

void pars_zero_one(int n, const double *p, fitpars *pars)
{
  int i;

  /* Translate parameters to real numbers in the model */
  assert(n <= pars->n);
  for (i=0; i < n; i++) {
    // Store the [0,1] values 
    pars->value[i] = p[i];
    // Store the real values
    *(pars->address[i]) = pars_from_01(pars,i,p[i]);
  }
}

void pars_set_constrained(fitpars *pars)
{
  int i;

  /* Translate real numbers in the model to parameters */
  pars->mode = PARS_CONSTRAINED;
  for (i=0; i < pars->n; i++) {
    double *p = pars->address[i], min=pars->min[i];
    if (*p < min) pars->value[i] = min;
    else if (*p > min+pars->range[i]) pars->value[i] = min+pars->range[i];
    else pars->value[i] = *p;
  }
}

void pars_constrained(int n, const double *p, fitpars *pars)
{
  int i;

  /* Translate parameters to real numbers in the model */
  assert(n <= pars->n);
  for (i=0; i < n; i++) {
    double v = p[i], min=pars->min[i];
    if (v < min) v = min;
    else {
      double max=min+pars->range[i];
      if (v > max) v = max;
    }
    *(pars->address[i]) = v;
  }
}

void pars_set_range(fitpars *pars, int ipar, double newmin, double newmax) 
{
  // Set the new minimum
  pars->min[ipar] = newmin;
  // Set the new range
  pars->range[ipar] = newmax-newmin;
  // Set the value between zero and one
  if (*(pars->address[ipar]) < pars->min[ipar]) pars->value[ipar] = 0.;
  else if (*(pars->address[ipar]) > pars->min[ipar]+pars->range[ipar]) pars->value[ipar] = 1.;
  else pars->value[ipar] = (*(pars->address[ipar])-pars->min[ipar])/pars->range[ipar];
}

/* $Id$ */
