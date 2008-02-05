/* This program is public domain. */

#ifndef _PARS_H
#define _PARS_H

/* FIXME this is what we should do rather than what we are
 * currently doing.
 */

/* The pars structure is designed to hold pointers to values that
 * can be varied by a fit.  Usually these will be held within a
 * model structure, but they can also be global control parameters,
 * etc.  This allows the fit program to take a simple vector of
 * doubles and not have to deal with issues such as which parameters
 * are active or where they are stored.
 *
 * pars_init(pars)
 *   prep the pars structure
 * pars_reset(pars)
 *   clear any existing parameters
 * pars_destroy(pars)
 *   clear the pars structure and free any associated memory
 * pars_print(pars)
 *   display the list of parameters
 * pars_print_set(pars,set)
 *   display a saved set of parameters without inserting them in the model
 * pars_count(pars)
 *   return the number of parameters
 * pars_vector(pars)
 *   return a work vector big enough to hold all parameters; there
 *   is only one copy of the vector and no locking mechanism so
 *   keep access short.  pars_add, pars_destroy and pars_extend
 *   will change the pointer.
 *
 *
 * Each parameter has:
 *    name    - string name of the parameter
 *    address - location in model to store parameter values
 *    min,max - box constraints on parameter values
 *
 * i = pars_add(pars,name,address,min,max)
 *   add a new parameter, returning its number
 *
 * pars_name(pars,i)
 * pars_min(pars,i)
 * pars_max(pars,i)
 * pars_range(pars,i)
 *   return the info for parameter i, 0-origin
 *
 * pars_set_range(pars,i,min,max)
 * pars_set_address(pars,i,address)
 *   reset info form parameter i, 0-origin
 *
 * pars_get(pars,v)
 * pars_set(pars,v)
 *   get/set the value of the parameter from the model into a vector
 *   of doubles.
 * v = pars_peek(pars,i)
 * pars_poke(pars,i,v)
 *   get/set the value of a single parameter from the model.
 *   
 * v = pars_to_01(pars,i,v)
 * v = pars_from_01(pars,i,v)
 *   convert v to/from [0,1] range for parameter i, 0-origin.
 * pars_get01(pars,v)
 * pars_set01(pars,v)
 * v = pars_peek01(pars,i)
 * pars_poke01(pars,i,v)
 *   same, but using scale independent box constraints which map from 
 *   values from [0,1] in the search space onto the entire parameter
 *   range in the model space.
 *
 * i = pars_select(pars,str)
 *   command-line command to select a parameter.  If str is provided,
 *   try to get the value from there first, otherwise let the user
 *   enter a parameter number, and display a list of available parameters
 *   if needed.
 *
 * pars_extend(pars,n)
 *   make sure there is room for at least n more parameters.  pars_add
 *   does this itself, but if you know you have a whole lot of parameters,
 *   use pars_extend first to keep from thrashing memory.
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct fitpars_struct {
  int n, capacity;
  double **address;
  const char ** name;
  double *min, *range, *value;
  int mode;
} fitpars;

double pars_to_01(const fitpars *pars, int i, double v);
double pars_from_01(const fitpars *pars, int i, double v);
void pars_init(fitpars *pars);
void pars_reset(fitpars *pars);
void pars_destroy(fitpars *pars);
int pars_extend(fitpars *pars, int n);

/* These will go away */
void pars_set_zero_one(fitpars *pars);
void pars_zero_one(int n, const double *p, fitpars *pars);
void pars_set_constrained(fitpars *pars);
void pars_constrained(int n, const double *p, fitpars *pars);


void pars_add(fitpars *pars, const char *name, 
	      double *d, double min, double max);

double* pars_vector(fitpars *pars);
int pars_count(const fitpars *pars);
const char* pars_name(const fitpars *pars, int n);
double pars_min(const fitpars *pars, int i);
double pars_max(const fitpars *pars, int i);
double pars_peek(const fitpars *pars, int i);
double pars_peek01(const fitpars *pars, int i);
void pars_poke(const fitpars *pars, int i, double v);
void pars_poke01(const fitpars *pars, int i, double v);
void pars_set(const fitpars *pars, const double v[]);
void pars_set01(const fitpars *pars, const double v[]);
void pars_get(const fitpars *pars, double v[]);
void pars_get01(const fitpars *pars, double v[]);

void pars_print(fitpars *pars);
int pars_select(fitpars *pars);
int pars_enter_value(fitpars *pars, int ipar, double *v);
void pars_set_range(fitpars *pars, int ipar, double newmin, double newmax);

#define PARS_ZERO_ONE 0
#define PARS_CONSTRAINED 1

#ifdef __cplusplus
}
#endif

#endif /* _PARS_H */

/* $Id$ */
