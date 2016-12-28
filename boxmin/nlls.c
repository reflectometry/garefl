/** \file
 * The Levenberg-Marquardt fitting routines.
 */

#include <stdio.h>
#include <math.h>
#include "lm.h"
#include "nlls.h"

#ifdef USE_SINGLE
# define levmar_dif slevmar_dif
# define levmar_bc_dif slevmar_bc_dif
#else
# define levmar_dif dlevmar_dif
# define levmar_bc_dif dlevmar_bc_dif
#endif

// like nlls_fn without the const qualifiers
typedef void (*lmfn)(Real *p, Real *hx, int m, int n, void *data);


void nlls_printstop(Real info[])
{
  switch ((int)(rint(info[6]))) {
  case 1: printf("small gradient\n"); break;
  case 2: printf("small Dp\n"); break;
  case 3: printf("itmax reached\n"); break;
  case 4: printf("singular matrix; restart with increased mu\n"); break;
  case 5: printf("stuck; restart with increased mu\n"); break;
  case 6: printf("small ||e||\n"); break;
  default: printf("unknown stopping criterion: info[6]=%g\n",info[6]); break;
  }
}

int nlls_worksize(int n, int k)
{
  return LM_DIF_WORKSZ(n, k);
}

void nlls(nlls_fn f, int n, int k, 
	  void *data, const Real y[], Real p[], Real covar[])
{
  Real info[LM_INFO_SZ];
  Real opts[LM_OPTS_SZ];
  opts[0]=LM_INIT_MU; /* Initial scale */
  opts[1]=opts[2]=opts[3]=1e-4; /* Stopping threshhold */
  opts[4]=LM_DIFF_DELTA; /* dx for numerical differentiation */
  levmar_dif((lmfn)f,
	      p, (Real *)y, n, k,
	      500,  /* itmax */
	      opts, /* opts[5] = { mu, e1, e2, e3, delta } */
	      info, /* info[LM_INFO_SZ] */
	      NULL, /* work vector of size LM_DIF_WORKSZ(3,k) */
	      covar, data);
  nlls_printstop(info);
}

void box_nlls(nlls_fn f, int n, int k, 
	      void *data, const Real y[],
	      const Real bounds[], Real p[], Real covar[])
{
  Real info[LM_INFO_SZ];
  Real opts[LM_OPTS_SZ];
  opts[0]=LM_INIT_MU; /* Initial scale */
  opts[1]=opts[2]=opts[3]=1e-4; /* Stopping threshhold */
  opts[4]=LM_DIFF_DELTA; /* dx for numerical differentiation */
  levmar_bc_dif((lmfn)f,
		 p, (Real *)y, n, k,
		 (Real *)bounds, (Real *)bounds+n,
		 500,  /* itmax */
		 opts, /* opts[5] = { mu, e1, e2, e3, delta } */
		 info, /* info[LM_INFO_SZ] */
		 NULL, /* work vector of size LM_DIF_WORKSZ(3,k) */
		 covar, data);
  nlls_printstop(info);
}
