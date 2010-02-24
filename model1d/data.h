/* This program is public domain. */

#ifndef _DATA_H
#define _DATA_H


#ifdef __cplusplus
extern "C" {
#endif

/* Accumulate weighted differences between model and theory. */
void
wsumsq(const int n, const Real x[], const Real y[], const Real dy[],
       const int Tn, const Real Tx[], const Real Ty[], const Real Tdy[],
       int *df, Real *sumsq);

typedef struct fitdata_struct {
  int capacity;
  int n;
  int have_uncertainty, have_resolution;
  Real *Q, *dQ, *R, *dR;
  char *file;
  char *comment;
} fitdata;

  /* Error string associated with error code returned from data_load */
const char *data_error(int code);
  /* Print the data to the screen */
void data_print(const fitdata *data);
  /* Initialize data structure for data */
void data_init(fitdata *data);
  /* Free memory associated with data */
void data_destroy(fitdata *data);
  /* Make sure the data can hold at least n elements */
int data_create(fitdata *data, int n);
  /* Load the data file, returning code=0 if success or code<0 if failure.
   * Call data_error to translate the returned code into a string.
   */
int data_load(fitdata *data, const char *file);
  /* Convert data from log to linear if it is substantially negative.
   * Return true if data was converted. */
int data_log2lin_if_negative(fitdata *data);
  /* Convert data from log to linear. */
void data_log2lin(fitdata *data);
  /* Print Q dQ R dR and theory line to a named file; the theory
     line is assumed to be sampled at the points in data.
   */
int data_printfit(const char *file, const fitdata *data, const Real fitR[]);
  /* Add partial chisq and degrees of freedom information from the fit
   * to the data file.  The fit points may be a subset of the total
   * points in the data file.  The points are not weighted by the uncertainty
   * in the data.
   */
void data_sumsq(const fitdata *data,
		const int fitn, const Real fitQ[], const Real fitR[],
		int *n, Real *sumsq);
  /* Add partial chisq and degrees of freedom information from the fit
   * to the data file.  The fit points may be a subset of the total
   * points in the data file.  The points are weighted by the uncertainty
   * in the data.
   */
void data_wsumsq(const fitdata *data,
		 const int fitn, const Real fitQ[], const Real fitR[],
		 int *n, Real *sumsq);

  /* Find the data points at which the theory needs to be calculated.
   * If Q is present for any cross section, then R needs to be calculated
   * for all cross-sections, so this is a question of merging the Q values
   * for all cross-sections into one set.
   * data_countQ() returns the number of values that will result.
   * data_mergeQ() merges the cross sections.
   */
int data_countQ(const fitdata *A, const fitdata *B,
		const fitdata *C, const fitdata *D);
void data_mergeQ(const fitdata *A, const fitdata *B,
		 const fitdata *C, const fitdata *D, Real Q[]);

/* Assign resolution to points in the data.  Resolution is assigned
 * in Q ranges, with the same resolution for -Q as Q.  Use Qhi == 0
 * to go to the maximum range in Q.
 *
 * Fixed slits:
 *   data_resolution_fixed(data,L,dLoL,0.,0.,dT);
 *
 * Opening slits:
 *   data_resolution_varying(data,L,dLoL,0.,0.,dToT);
 *
 * Fixed below |Qlo| then opening:
 *   data_resolution_fv(data,L,dLoL,|Qlo|,dTlo,dToT);
 *
 * Opening between |Qlo| and |Qhi|, fixed above and below:
 *   data_resolution_fvf(data,L,dLoL,|Qlo|,|Qhi|,dTlo,dToT,dThi);
 *
 * Constant resolution (e.g., from TOF source)
 *   data_resolution(data, res);
 *
 * If instrument uses angle T rather than Q:
 *   #include "reflcalc.h"
 *   Q = T2Q(L,T)
 *
 * To compute dT from slit openings and separation use:
 *   #include "reflcalc.h"
 *   dT = resolution_dT(s1,s2,d)
 *
 * For opening slits, dToT also needs an incident angle in degrees:
 *   #include "reflcalc.h"
 *   dToT = resolution_dToT(s1,s2,d,T)
 *
 * For time of flight instruments with binning adjusted to maintain
 * constant dQoQ over a particular Q range:
 *   data_resolution_dQoQ(data, dQoQ, |Qlo|, |Qhi|) // Part range
 *   data_resolution_dQoQ(data, dQoQ, 0., 0.)       // Full range
 */
void
data_resolution_fixed(fitdata *data, Real L, Real dLoL,
		      Real Qlo, Real Qhi, Real dT);
void
data_resolution_varying(fitdata *data, Real L, Real dLoL,
			Real Qlo, Real Qhi, Real dToT);

void
data_resolution_fv(fitdata *data, Real L, Real dLoL,
		   Real Qlo, Real dTlo, Real dToT);

void
data_resolution_fvf(fitdata *data, Real L, Real dLoL,
		    Real Qlo, Real Qhi,
		    Real dTlo, Real dToT, Real dThi);
void
data_resolution_dQoQ(fitdata *data, Real dQoQ, Real Qlo, Real Qhi);

void
data_constant_resolution(fitdata *data, Real res);

#ifdef __cplusplus
}
#endif

#endif

/* $Id$ */
