/* This program is public domain. */

#ifndef _AMOEBA_H
#define _AMOEBA_H
/** \file
   Amoeba bounds-constrained nelder-mead simplex.

   Use the following to fit function f(np,p,data) starting from p=Po

	 \code
   // space for the simplex
   simplex s;
   p = malloc(sizeof(double)*(n*n+4*n+1));

   // initialize the simplex
   amoeba_init(&s,n,bounds,p,f,data);
   amoeba_reset(&s,Po,10.);

   // do the fit
   pbest = amoeba(&s,tolerance,max_iterations);

   // print the result
   printf("minimum %g found at (%g",pbest[n],pbest[0]);
   for (i=1; i < n; i++) printf(",%g",pbest[i]);
   printf(")\n");
	 \endcode

You can use your own outer loop as well:

   \code
   pbest = amoeba_best(&s);
   for (it=0; i < max_iterations; it++) {
     // Determine if the simplex is stuck in a flat region
     if (amoeba_flatness(&s) < tolerance) break;
     pbest = amoeba_step(&s);
   }
	 \endcode

Use the following macros to access the simplex directly:

   \code
   double *pk = amoeba_VERTEX(&s,k);
   amoeba_VALUE(&s,k) = fn(n,pk,userdata);
	 \endcode

You may want to do this if, for example, you were
filling the simplex from the current population of
a genetic algorithm for which you have already computed
the fitness of each vertex.

*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

	/** Objective function typedef.  Your objective functions should be of
			this type to minimize them with boxmin.
	*/
typedef double (*optimfn)(int, const double[], void *);
typedef struct simplex__ 
{
  int n;          /* n dimensions */
  double *p;      /* n+1 vertices of the simplex, n coordinates for each */
  double *ptry;   /* a trial vertex with n coordinates */
  double *psum;   /* sum of the vertex coordinates across all vertices */
  const double *bounds; /* NULL or 2*n bounds (lo 1 2 3 ... hi 1 2 3 ...) */
  int ilo, ihi, inhi; /* lowest and highest vertices */
  int iterations; /* Number of iterations */
  optimfn fn;     /* Function to optimize */
  void *userdata; /* Extra data for objective function */
} simplex;

	/** returns the number of doubles needed: n^2+4n+1 */
int amoeba_worksize(int n);

	/** set up the simplex structure for the problem. */
void amoeba_init(
 simplex *s,            /* structure to hold the state of the simplex. */
 int n,                 /* the number of fit variables. */
 const double bounds[], /* the lo-hi bounds for each variable. */
 double p[],            /* the simplex workspace (min. n^2+4n+1 values) */ 
 optimfn fn,          /* f(int np, double p[], void *data) */
 void *userdata);     /* extra data for function f */

	/** Some macros for accessing vertices and values directly in the simplex */
#define amoeba_VERTEX(s,k) ((s)->p+((s)->n+1)*(k))
#define amoeba_VALUE(s,k) ((s)->p)[((s)->n+1)*(k+1)-1]

	/** set initial vertices from Po using the metropolis structure
 * If amoeba is bound-constrained, scale is a percentage [0 to 100]
 * of the initial search space to jump in each direction.
 *
 * \param s Simplex returned by call to amoeba()
 * \param Po Initial value for the minimization
 * \param scale Size of the initial Simplex (as portion of search space).
 * \returns Nothing
 */
void amoeba_reset(simplex *s, double Po[], double scale);

	/** print the current simplex to the screen */
void amoeba_dumpsimplex(simplex *s);

	/** Search the vertices, remembering and returning the best */
double* amoeba_best(simplex *s);

	/** Take one amoeba step, recalculating any vertices as necessary.
 * Returns the new best vertex.
 */
double* amoeba_step(simplex *s);

	/** return the flatness of the simplex
 * Call this after amoeba or amoeba_step.
 */
double amoeba_flatness(simplex *s);

double* amoeba(simplex *s, double ftol, int itmax);

#ifdef __cplusplus
}
#endif

#endif /* _AMOEBA_H */

/* $Id: amoeba.h 35 2007-05-25 15:12:45Z ziwen $ */
