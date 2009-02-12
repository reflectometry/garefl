#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

/** \file
 * Compile with -DCLAMP_ENDS=0 to reproduce pgs/matlab behaviour driving
 * the splines to 0 beyond the ends.  The default CLAMP_ENDS=1 sets the
 * value to the value of the final knot.
 */

#ifndef CLAMP_ENDS
#define CLAMP_ENDS 1
#endif

#ifdef HAVE_INLINE
inline double min(double a, double b) { return a < b ? a : b; }
inline double max(double a, double b) { return a > b ? a : b; }
#else
#define min(a,b) ( (a) < (b) ? (a) : (b) )
#define max(a,b) ( (a) > (b) ? (a) : (b) )
#endif

#include<stdio.h>
#include "reflconfig.h"

/* function prototype */
static int
find_interval(const int n, const double x[], const double v);

/* trisolve(n, l[],d[],u[], b[])
 * Solves the tridiagonal system A*x = b.  The vector l is the subdiagonal,
 * d is the diagonal and u is the superdiagonal of A.  On return the vector
 * B contains the solution.
 *
 * l,d,u are modified during the calculation.
 *
 * Uses LAPACK's nonsymmetric tridiagonal solver dgtsv.
 */
#define Fdgtsv F77_FUNC(dgtsv,DGTSV)
extern int Fdgtsv(const int *n, const int *nrhs, const double l[],
		  const double d[], const double u[], double b[],
		  const int *bskip, int *info);
int trisolve(const int n, double l[], double d[], double u[], double b[])
{
  int info, bskip=n, nrhs=1;
  Fdgtsv(&n, &nrhs, l, d, u, b, &bskip, &info);
  return info == 0;
}

/* Return the control value for a segment.  Beyond the range of control
 * values, this is either clamped to the last control value if CLAMP_ENDS==1
 * or to 0 if CLAMP_ENDS==0.
 */
inline double
control_value(const int n, const double c[], const int k)
{
#if CLAMP_ENDS
  return k < 0 ? c[0] : (k >= n ? c[n-1] : c[k]);
#else
  return k < 0  || k >= n ? 0. : c[k];
#endif
}

/* Refine a B-spline, assuming it is uniform */
/* This uses the algorithm described in http://graphics.idav.ucdavis.edu/education/CAGDNotes/Cubic-B-Spline-Curve-Refinement/Cubic-B-Spline-Curve-Refinement.html (with the endpoints of the original control polygon duplicated) */
/* Given a sequence of n control points, will return a sequence of 2n+1 control points */

void
bspline3_refine(int n, const double control[], double new_points[])
{
  /*printf("N is %d\n", n);*/
  int i, j, k;
  double edge_points[n+1];
  double vertex_points[n];
  edge_points[0] = control[0];
  for(i=1; i < n; i++)
    edge_points[i] = ((control[i-1] + control[i])/2.0);
  edge_points[n] = control[n-1];
  for(j=0; j < n; j++)
    vertex_points[j] = (edge_points[j] + (2.0 * control[j]) + edge_points[j+1]) / 4.0;
  new_points[0] = edge_points[0];
  for(k=0; k < n; k++)
    {
      new_points[(2*k)+1] = vertex_points[k];
      new_points[(2*k)+2] = edge_points[k+1];
    }
}

/* Evaluate a B-spline at a set of points within a segment.  */

static void
eval_points(const int N_knots,
	    const double knot[], const double control[],
	    const int segment, int n, const double t[],
	    double v[])
{
  int N_control = N_knots - 4;
  int i;
  const double tm2 = knot[max(segment-2,0)];
  const double tm1 = knot[max(segment-1,0)];
  const double tm0 = knot[segment-0];
  const double tp1 = knot[min(segment+1,N_knots-1)];
  const double tp2 = knot[min(segment+2,N_knots-1)];
  const double tp3 = knot[min(segment+3,N_knots-1)];

  const double P4o = control_value(N_control, control, segment-0);
  const double P3o = control_value(N_control, control, segment-1);
  const double P2o = control_value(N_control, control, segment-2);
  const double P1o = control_value(N_control, control, segment-3);

  /* printf("P1 %g P2 %g P3 %g P4 %g\n",P1,P2,P3,P4); */
  for (i=0; i < n; i++) {
    double ti = t[i];

    double P4 = P4o;
    double P3 = P3o;
    double P2 = P2o;
    double P1 = P1o;

    P4 = ( (ti-tm0)*P4 + (tp3-ti)*P3 ) / (tp3 - tm0);
    P3 = ( (ti-tm1)*P3 + (tp2-ti)*P2 ) / (tp2 - tm1);
    P2 = ( (ti-tm2)*P2 + (tp1-ti)*P1 ) / (tp1 - tm2);

    /* printf("ti %g: P2 %g P3 %g P4 %g\n",ti,P2,P3,P4); */


    P4 = ( (ti-tm0)*P4 + (tp2-ti)*P3 ) / (tp2 - tm0);
    P3 = ( (ti-tm1)*P3 + (tp1-ti)*P2 ) / (tp1 - tm1);

    P4 = ( (ti-tm0)*P4 + (tp1-ti)*P3 ) / (tp1 - tm0);

    v[i] = P4;
  }
}


/* Evaluate a B-spline at a single point within a segment. */
inline double
eval_point(const int N_knots,
	   const double knot[], const double control[],
	   const int segment, const double ti)
{
  const double tm2 = knot[max(segment-2,0)];
  const double tm1 = knot[max(segment-1,0)];
  const double tm0 = knot[segment-0];
  const double tp1 = knot[min(segment+1,N_knots-1)];
  const double tp2 = knot[min(segment+2,N_knots-1)];
  const double tp3 = knot[min(segment+3,N_knots-1)];
  const int N_control = N_knots - 4;

  double P4 = control_value(N_control, control, segment-0);
  double P3 = control_value(N_control, control, segment-1);
  double P2 = control_value(N_control, control, segment-2);
  double P1 = control_value(N_control, control, segment-3);

  P4 = ( (ti-tm0)*P4 + (tp3-ti)*P3 ) / (tp3 - tm0);
  P3 = ( (ti-tm1)*P3 + (tp2-ti)*P2 ) / (tp2 - tm1);
  P2 = ( (ti-tm2)*P2 + (tp1-ti)*P1 ) / (tp1 - tm2);

  P4 = ( (ti-tm0)*P4 + (tp2-ti)*P3 ) / (tp2 - tm0);
  P3 = ( (ti-tm1)*P3 + (tp1-ti)*P2 ) / (tp1 - tm1);

  P4 = ( (ti-tm0)*P4 + (tp1-ti)*P3 ) / (tp1 - tm0);

  return P4;
}


/* Evaluate a B-spline and all derivatives at a point within a segment. */
void
bspline3_eval_all_derivs(const int N_knots,
		const double knot[], const double control[],
		const double t, double derivs[4])
{
  const int N_control = N_knots - 4;
  const int segment = find_interval(N_knots, knot, t);

  if (segment >= N_knots) {

    derivs[0] = control_value(N_control, control, segment);
    derivs[1] = derivs[2] = derivs[3] = 0;
#if 0
    printf("t=%g tmax=%g Cmax=%g ret=%g\n",
	   t,knot[N_knots-1],control[N_knots-5],derivs[0]);
#endif

  } else {

    const double tm2 = knot[max(segment-2,0)];
    const double tm1 = knot[max(segment-1,0)];
    const double tm0 = knot[segment-0];
    const double tp1 = knot[min(segment+1,N_knots-1)];
    const double tp2 = knot[min(segment+2,N_knots-1)];
    const double tp3 = knot[min(segment+3,N_knots-1)];

    double P4 = control_value(N_control, control, segment-0);
    double P3 = control_value(N_control, control, segment-1);
    double P2 = control_value(N_control, control, segment-2);
    double P1 = control_value(N_control, control, segment-3);

    double Q4 = (P4 - P3) * 3 / (tp3 - tm0);
    double Q3 = (P3 - P2) * 3 / (tp2 - tm1);
    double Q2 = (P2 - P1) * 3 / (tp1 - tm2);

    double R4 = (Q4 - Q3) * 2 / (tp2 - tm0);
    double R3 = (Q3 - Q2) * 2 / (tp1 - tm1);

    double S4 = (R4 - R3) * 1 / (tp1 - tm0);

#if 0
    printf("segment=%d, max=%d\n",segment,N_knots);
    printf("t: %f %f %f %f %f %f\n", tm2, tm1, tm0, tp1, tp2, tp3);
    printf("P: %f %f %f %f\n", P4,P3,P2,P1);
    printf("Q: %f %f %f\n", Q4,Q3,Q2);
    printf("R: %f %f\n", R4,R3);
    printf("S: %f\n", S4);
#endif

    P4 = ( (t-tm0)*P4 + (tp3-t)*P3 ) / (tp3 - tm0);
    P3 = ( (t-tm1)*P3 + (tp2-t)*P2 ) / (tp2 - tm1);
    P2 = ( (t-tm2)*P2 + (tp1-t)*P1 ) / (tp1 - tm2);

    P4 = ( (t-tm0)*P4 + (tp2-t)*P3 ) / (tp2 - tm0);
    P3 = ( (t-tm1)*P3 + (tp1-t)*P2 ) / (tp1 - tm1);

    P4 = ( (t-tm0)*P4 + (tp1-t)*P3 ) / (tp1 - tm0);

    Q4 = ( (t-tm0)*Q4 + (tp2-t)*Q3 ) / (tp2 - tm0);
    Q3 = ( (t-tm1)*Q3 + (tp1-t)*Q2 ) / (tp1 - tm1);

    Q4 = ( (t-tm0)*Q4 + (tp1-t)*Q3 ) / (tp1 - tm0);

    R4 = ( (t-tm0)*R4 + (tp1-t)*R3 ) / (tp1 - tm0);

    derivs[0] = P4;
    derivs[1] = Q4;
    derivs[2] = R4;
    derivs[3] = S4;
  }

  /* FIXME If performance is an issue try defining constants for t-tk,
   * and building a table of constants for each segment i:
   *    P_i,Q_i,R_i,S_i, t_{i+1}-t_i, t_{i+2}-t_i and t_{i+3}-t_i
   * This will save a number of subtractions and a few divisions per node.
   * Set up the tables and the calls so that segment lies in [3,knots-3]
   * and no tests are required. */
}


/* Evaluate a number of grid points within a B-Spline segment.
 * The evaluation uses the t_values given in the argument.
 * It is bounded by the edge of the segment.
 *
 * The segment should be in the range [-1,N_knots-1], inclusive.
 *
 * Returns the (t,v) values of the grid points and the number of
 * points evaluated.
 */
static int
eval_grid(const int N_knots,
	  const double knot[], const double control[],
	  const int segment, const double t_values[],
	  const int n, double t[], double v[])
{
  int i,k;
  assert(segment>=-1 && segment<N_knots);

  /* Determine the number of points to evaluate. */
  if (segment >= N_knots) {
    k = n; /* After the spline, evaluate all */
  } else if (t_values[0] > knot[segment+1]) {
    return 0;
  } else {
    for(k=0; t_values[k+1] < knot[segment+1]; k++)
      ;
    k++;
    if (k > n) k = n;
  }

  /* Set up evaluation points */
  /* FIXME restore 'grid' behaviour of this function and use a
   * separate function for evaluating a set of non-uniform points
   * in a segment.
   */
  for (i=0; i < k; i++) t[i] = t_values[i];

  /* Handle points beyond the ends */
  if (segment < 0 || segment >= N_knots) {
    double vk = control_value(N_knots-4,control,segment);
    for (i=0; i < k; i++) v[i] = vk;
    return k;
  }

#if 1

  /* Evaluate all the points separately.  Do this when accuracy is
     more important than speed. */
  if (k) eval_points(N_knots, knot, control, segment, k, t, v);

#else

  /* If only need a few points, evaluate them directly */
  if (k < 4) {
    if (k) eval_points(N_knots, knot, control, segment, k, t, v);
    return k;
  } else {
    double d[4];
    double z,dz1,dz2,dz3;

    eval_all_derivs(N_knots, knot, control, segment, t[0], d);
    z = d[0]; dz = d[1]; d2z = d[2]; d3z= d[3];

    /* Step along the curve using the current derivative.  Update the
     * derivative using the second derivative.  Update the second derivative
     * using the third derivative.  This is a cubic, so the third derivative
     * is a constant and doesn't need to be updated. */
    for (i=1; i < k; i++) {
      d2z += d3z;
      dz += d2z;
      z += dz;
      v[i] = z;
    }
  }
#endif

  return k;
}

/* Determine which interval of x contains the value v.
 * x is assumed to be ascending and have more than one element.
 * Returns the interval containing the value, -1 if before the
 * first or n if after the last.  For blank intervals, return
 * the last in the series.
 */
static int
find_interval(const int n, const double x[], const double v)
{
  int lo = 0, hi=n-1;
  if (v < x[0]) return -1;
  if (v >= x[n-1]) return n;
  while (lo < hi-1) {
    const int mid = (lo+hi)/2;
    /* Change '<' to '<=' below if you want to return the first in the
     * a series of identical knots rather than the last. You will
     * also need to change the '<' limits test above. */
    if (v < x[mid]) hi=mid;
    else lo=mid;
  }
  return lo;
}

/* Fill a grid of values.  Complete docs in bspline.h. */
void
bspline3(const int N_knots,
	 const double knot[], const double control[],
	 const double t_vals[],
	 const int n, double t[], double v[])
{
  int k, segment;

  /* Process segments until the entire mesh is filled. */
  segment = find_interval(N_knots, knot, t_vals[0]);
  k = eval_grid(N_knots, knot, control, segment, t_vals, n, t, v);
  while (k < n) {
    k += eval_grid(N_knots, knot, control, ++segment,
		   t_vals+k, n-k, t+k, v+k);
    /* The loop is guaranteed to halt because eval_grid fills to the end
     * beyond the last knot.  But a little paranoia never hurt anybody...
     */
    assert(segment < N_knots+1);
  }
}

/* Evaluates the derivatives of a spline defined be N_knots, knot, and
 * control at the given currT and returns the result in
 * result[]. result[] is in the format [a, b, c, d] where a = f'(t), b
 * = f''(t), c = f'''(t), and d = the distance from t to the next
 * knot. Note that this function may call eval_points with a t-value
 * that is in the "wrong" segment (e.g. calling eval_points with
 * segment = 0 but where the T-value is not in fact in the first
 * segment) and eval_points needs to work properly (i.e. return the
 * value of the function describing that segment for the given
 * t-value, which will not in general be a point actually on the
 * spline) in those cases. */

void
bspline3_eval_derivs(const int N_knots, const double knot[],
	    const double control[], const double currT,
		     double result[])
{
  int segment = 0;
  double values[4];
  double tvals[4] = {currT, currT+1.0, currT+2.0, currT+3.0};
  segment = find_interval(N_knots, knot, currT);
  eval_points(N_knots, knot, control, segment, 4, tvals, values);
  result[0] = (-1.833333*values[0] + 3.0*values[1] - 1.5*values[2] + 0.333333*values[3]);
  result[1] = (2.0*values[0] - 5.0*values[1] + 4.0*values[2] - 1.0*values[3]);
  result[2] = (-1*values[0] + 3.0*values[1] - 3.0*values[2] + 1.0*values[3]);
  result[3] = knot[segment+1] - currT;
}

/* Interpolates a spline through the data points.
 * Takes N_data data points and returns a set of N_data+2
 * control points that create a spline passing through the data points.
 * Requires a work vector of length (3*N_data)+4.
 */

int
bspline3_interpolate(const int N_data, double data[],
		     double control[], double work[])
{
  int i;

  /* FIXME need test cases. */
  double* bvect = control; /*N_data + 2*/
  double* ldiag = work; /*N_data + 1*/
  double* mdiag = work + (N_data + 1); /*N_data + 2*/
  double* udiag = work + (2*N_data + 3); /*N_data + 1*/

  ldiag[0] = -0.3;
  ldiag[1] = 0.25;
  for(i=2; i<(N_data - 1); i++)
    ldiag[i] = (1.0 / 6.0);
  ldiag[N_data - 1] = 0;
  ldiag[N_data] = 0;

  udiag[0] = 0.0;
  udiag[1] = 0.0;
  for(i=2; i<(N_data - 1); i++)
    udiag[i] = (1.0 / 6.0);
  udiag[N_data - 1] = 0.25;
  udiag[N_data] = 0.3;

  mdiag[0] = 1.0;
  mdiag[1] = 0.3;
  mdiag[2] = (7.0 / 12.0);
  for(i=3; i<(N_data-1); i++)
    mdiag[i] = (2.0 / 3.0);
  mdiag[N_data - 1] = (7.0 / 12.0);
  mdiag[N_data] = -0.3;
  mdiag[N_data + 1] = 1.0;

  bvect[0] = data[0];
  /* bvect[1] = (data[1] - data[0]) * N_data; */ /* estimate dx/dt at beginning */
  bvect[1] = 0;
  memcpy(&bvect[2], &data[1], (N_data - 2)*(sizeof(double)));
  /* bvect[N_data] = (data[N_data - 1] - data[N_data - 2]) * N_data; */ /* estimate dx/dt at end */
  bvect[N_data] = 0;
  bvect[N_data + 1] = data[N_data - 1];

  /* Note: trisolve destroys the values during the computation. */
  /* Returns false if the system is not solvable. */
  if (!trisolve(N_data + 2, ldiag, mdiag, udiag, bvect)) return 0;

  return 1;
}


#ifdef TEST

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Print an error if a is not within tolerance of b. */
void check(double a, double b, double tol)
{
  if (fabs(a-b) > tol) {
    fprintf(stderr,"%g != %g within tolerance %g [a-b=%g]\n",a,b,tol,fabs(a-b));
  }
}

/* Print an error if vector a is not within tolerance of b. */
void checkv(int k, const double a[], const double b[], double tol)
{
  int i, error=0;
  for (i=0; i < k; i++) error += (fabs(a[i]-b[i]) > tol);
  if (error) {
    printf("vectors differ in %d values:\n",error);
    for (i=0; i < k; i++)
      printf("%d: %c %g %g %g\n", i, fabs(a[i]-b[i])>tol?'*':' ', a[i], b[i],
		      b[i]-a[i]);
  }
}

#define NT 10100
int main(int argc, char *argv[])
{
  const double knot1[] = { 0, 0, 0, 1, 1, 3, 4, 6, 6, 6 };
  const double knot2[] = { 0, 0, 0, 0, 1, 4, 5, 5, 5, 5 };
  const double knot3[] = { -1, 0, 0, 1, 1, 3, 4, 6, 6, 7 };
  const double c1[] = { 0, 0, 0, 0, 0, 0 };
  const double c2[] = { 1, 1, 1, 1, 1, 1 };
  const double c3[] = { 1, 2, 3, 4, 5, 6 };
  const double c4[] = {1e-5, 2e-5, 3e-5, 2e-5, 1e-5, 2e-5};
  const double t0 = knot3[0];
  double t[NT], v[NT];

  if (argc > 1) {
    int n = atoi(argv[1]);
    double dt = (knot3[9]-knot3[0])/n;
    int sumk = 0;
    while (sumk < n) {
      const int k = (n-sumk < NT ? n-sumk : NT);
      int i;
      bspline3(10, knot3, c4, t0+sumk*dt, dt, k, t, v);
      for (i=0; i < k; i++) printf("%.15g %.15g\n",t[i],v[i]);
      sumk += k;
    }
  } else {
    int i;
    double expected_t[] = {
      -1, -0.2, 0.6, 1.4, 2.2, 3, 3.8, 4.6, 5.4, 6.2, 7
    } ;
    double expected_v[] = {
      1, 1, 1.72, 2.58542222222222, 2.41973333333333, 1.84444444444444,
      1.41137777777778, 1.49257777777778, 1.8572, 2, 2
    } ;
    double dt;

    printf("# Checking single points\n");
    bspline3(10, knot1, c1, 2.2, 1., 1, t, v);
    check(v[0],0.,3e-16);
    bspline3(10, knot1, c2, 2.2, 1., 1, t, v);
    check(v[0],1.,3e-16);
    bspline3(10, knot2, c3, 2., 1., 1, t, v);
    check(v[0],761./240.,3e-15);
    bspline3(10, knot1, c3, 3.2, 1., 1, t, v);
    check(v[0],4.2976,3e-16);

    printf("# Checking sparse vector\n");
    dt = (knot3[9]-knot3[0])/10;
    bspline3(10, knot3, c4, t0, dt, 11, t, v);
    checkv(11,t,expected_t,5e-15);
    checkv(11,v,expected_v,5e-15);

    printf("# Checking single points\n");
    for (i=0; i < 11; i++) {
      bspline3(10, knot3, c4, expected_t[i], 0.001, 1, t, v);
      check(v[0],expected_v[i],5e-15);
    }

    dt = (knot3[9]-knot3[0])/100;
    bspline3(10, knot3, c4, t0, dt, 100, t, v);
    printf("# Checking dense vector t values\n");
    for (i=0; i < 11; i++) check(t[10*i],expected_t[i],5e-15);
    printf("# Checking dense vector v values\n");
    for (i=0; i < 11; i++) check(v[10*i],expected_v[i],5e-15);
  }

  return 0;
}

#endif
