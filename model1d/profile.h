/* This program is public domain. */

#ifndef _PROFILE_H
#define _PROFILE_H

#include "reflconfig.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Profiles are an array of simple rectangular slices.  They are produced
 * by a model and used to calculate ideal reflectivity.
 *
 * The following fields are part of a profile structure:
 *
 *    int length            - number of slices
 *    double vacuum_offset  - depth before the first interface
 *    double *d             - depth of each slice in Angstroms
 *    double *rho           - scattering length density of each slice in Nb
 *    double *mu            - absorption coefficient of each slice
 * #ifdef HAVE_MAGNETIC
 *    double *P             - magnetic scattering in Nb
 *    double *theta         - magnetic direction in degrees
 *    double *expth         - magnetic direction in rectangular coordinates
 * #endif
 *
 * Other fields may be present and the order may be different.
 * 
 * profile_init(profile*)
 *   Initialize a new profile structure.
 *
 * profile_destroy(profile*)
 *   Free memory associated with a profile structure.  This does not
 *   delete the structure.
 *
 * profile_reset(profile*)
 *   Empty the profile (set the length to 0) but do not free the
 *   memory associated with it.  Reset profile to a 'good' state.
 *   Use profile_destroy if you want to free memory.
 *
 * int profile_extend(profile*,int n)
 *   Extend the profile so that it can hold at least n more slices.
 *   Returns false if profile could not be extended.
 *
 * int profile_good(profile*)
 *   Return true if the profile is good, or false if there was a
 *   memory allocation error.  profile_reset resets the 'good' state.
 *
 * double profile_depth(profile*)
 *   Return the total depth of all slices in the profile, including
 *   slices in the vacuum interface.
 *
 * profile_slice(profile*,double d, double rho, double mu)
 *   Add a slice to the end of the profile. If HAVE_MAGNETIC, there
 *   are an additional two paramters: P and theta.
 *
 * profile_copy(profile*, int start, int length)
 *   Copy slices from start for length to the end of the profile.
 *
 * profile_expth(profile*)
 *   Convert theta from polar to rectangular coordinates.  The expth
 *   field is not maintained during profile construction.  Call 
 *   profile_expth after the profile is built.  Then you can use 
 *   p->expth as the expth argument to the magnetic reflectivity calculation.
 *
 * Profiles are built with dynamic memory.  Sometimes a new profile
 * will require more memory than is available, and in this case,
 * the profile is destroyed and an empty profile is returned.  You
 * can either test the build_profile return value, or you can check
 * if the profile it created is bad.
 *
 * Movie example:
 * 
 *   profile p;
 *   profile_init(&p);
 *   while (changes to model) {
 *      profile_reset(&p);
 *      model_profile(model*,&p);
 *      if (profile_good(&p)) {
 *         clear graph
 *         staircase(p.length,cumsum(p.d)-p.vacuum_offset,p.rho);
 *         staircase(p.length,cumsum(p.d)-p.vacuum_offset,p.mu);
 *      }
 *   }
 *   profile_destroy(&p);
 */

#ifdef HAVE_MAGNETIC
/* There are only 5 fields, but expth requires an extra two values */
#define PROFILE_FIELDS 7
#else
#define PROFILE_FIELDS 3
#endif

typedef struct profile_struct {
  int capacity; /* private---don't use it directly */
  int n;
  double vacuum_offset;
  double *rho;
  double *mu;
  double *d;
#ifdef HAVE_MAGNETIC
  double *P;
  double *theta;
  double *expth;
#endif
} profile;

void profile_init(profile *p);
void profile_reset(profile *p);
void profile_destroy(profile *p);
int profile_extend(profile *p, int n);
int profile_good(profile *p);
double profile_depth(profile *p);
void profile_copy(profile *p, int start, int length);
void profile_expth(profile *p);
void profile_slice(profile *p, double d, double rho, double mu
#ifdef HAVE_MAGNETIC
		   , double P, double theta
#endif
		   );
void profile_print(profile *p, char *file);

#ifdef __cplusplus
}
#endif

#endif /* _PROFILE_H */

/* $Id$ */
