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
        /* shift name pointer array */
	memmove(pars->address+need, pars->address+size, size*sizeof(void*));
      }
      if (pars->min != NULL) {
        /* shift value and range arrays */
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
    pars->value = pars->min + 2*need;
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
  return pars->name[i] ? pars->name[i] : "unknown";
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
  const double range = pars->range[i];
  const double min = pars->min[i];
  const double max = min + range;
  v = v < min ? min : (v > max ? max : v);
  *(pars->address[i])  = v;
}
void pars_poke01(const fitpars *pars, int i, double v)
{
  pars_poke(pars, i, pars_from_01(pars, i, v));
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


static void _print_one(fitpars *pars, int i, double v, int zero_one)
{
  char range[]="..........";
  double value, portion;
  if (zero_one) {
    portion = v;
    value = v*pars->range[i]+pars->min[i];
  } else {
    portion = (v-pars->min[i])/pars->range[i];
    value = v;
  }
  if (portion < 0.) portion = 0.; 
  else if (portion >= 1.) portion = 0.999999;

  printf("%3d ",i);
  if (pars->name[i] != NULL) printf("%25s", pars->name[i]);
  range[(int)floor(portion*strlen(range))] = '|';
  printf(" %s ", range);
  printf("%g in [%g,%g]\n",value,pars->min[i],pars->min[i]+pars->range[i]);
}

void pars_print_set(fitpars *pars, double *set, int zero_one)
{
  int i;
  for (i=0; i<pars->n; i++) {
    _print_one(pars, i, set[i], zero_one);
  }
}

void pars_print(fitpars *pars)
{
  int i;
  for (i=0; i<pars->n; i++) {
    _print_one(pars, i, pars_peek(pars,i), 0);
  }
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
    const char *name = pars_name(pars, p);
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

void pars_set_range(fitpars *pars, int ipar, double newmin, double newmax) 
{
  // Set the new minimum
  pars->min[ipar] = newmin;
  // Set the new range
  pars->range[ipar] = newmax-newmin;
  // Force the current value into the range.
  pars_poke(pars, ipar, pars_peek(pars, ipar));
}

/* $Id$ */
