
/** \file
 * The code to calculate the uncertainty ( uct ) to  a least squares fit
 */

/**===================================================================
 *  This is a code to calculate the uncertainty ( uct ) to  a least squares fit
 *
 *  All the notations and conventions are adopted from the work of
 *  Manolis Lourakis (lourakis@ics.forth.gr)
 *  Institute of Computer Science, Foundation for Research & Technology-Hellas
 *
 *   Ziwen Fu
 *   08/20/2007
 *========================================================================
 */

#ifndef LM_REAL 
#error This file should not be compiled directly!
#endif


/* precision-specific definitions */
#define LEVMAR_DIF_UN        LM_ADD_PREFIX(levmar_dif_un)
#define FDIF_FORW_JAC_APPROX LM_ADD_PREFIX(fdif_forw_jac_approx)
#define FDIF_CENT_JAC_APPROX LM_ADD_PREFIX(fdif_cent_jac_approx)
#define TRANS_MAT_MAT_MULT   LM_ADD_PREFIX(trans_mat_mat_mult)
#define AX_EQ_B_LU           LM_ADD_PREFIX(Ax_eq_b_LU_noLapack)




/* ===================================================================
 * This function computes the inverse of A ( m x m ). 
 * B = inv( A )
 *
 * Returns 0 in case of error, 1 if successful
 */
static int LM_LUINVERSE(LM_REAL *A,   /* INPUT  */
												LM_REAL *B,   /* OUTPUT */
												int m) {
  void *buf=NULL;

  register int i, j, k, l;
  int *idx, maxi=-1, idx_sz, a_sz, x_sz, work_sz, tot_sz;
  LM_REAL *a, *x, *work, max, sum, tmp;

  /* calculate required memory size */
  idx_sz = m;
  a_sz   = m*m;
  x_sz   = m;
  work_sz=m;
  tot_sz =idx_sz*sizeof(int) + (a_sz+x_sz+work_sz)*sizeof(LM_REAL);

  buf = (void *)malloc(tot_sz);
  if(!buf){
    fprintf(stderr, "memory allocation\n");
    exit(1);
  }

  idx  = (int *)buf;
  a    = (LM_REAL *)(idx + idx_sz);
  x    = a + a_sz;
  work = x + x_sz;

  /* copying A to a */
  for(i=0; i<a_sz; ++i)
		 a[i]=A[i];

  /* compute the LU decomposition of a row permutation of matrix a;
		 the permutation itself is saved in idx[] */
	for(i=0; i<m; ++i){
		  max=0.0;
		  for(j=0; j<m; ++j)
			   if((tmp=FABS(a[i*m+j]))>max)  max=tmp;
		     if(max==0.0){
         fprintf(stderr,  "Singular matrix\n");
         free(buf);
         return 0;
      }
		  work[i]=CNST(1.0)/max;
	}

	for(j=0; j<m; ++j){
		 for(i=0; i<j; ++i){
			  sum=a[i*m+j];
			  for(k=0; k<i; ++k)
           sum-=a[i*m+k]*a[k*m+j];
			  a[i*m+j]=sum;
		 }
		 max=0.0;
		 for(i=j; i<m; ++i){
			  sum=a[i*m+j];
			  for(k=0; k<j; ++k)
           sum-=a[i*m+k]*a[k*m+j];
			  a[i*m+j]=sum;
			  if((tmp=work[i]*FABS(sum))>=max){
				   max=tmp;
				   maxi=i;
			  }
		}
		if(j!=maxi){
			for(k=0; k<m; ++k){
				tmp=a[maxi*m+k];
				a[maxi*m+k]=a[j*m+k];
				a[j*m+k]=tmp;
			}
			work[maxi]=work[j];
		}
		idx[j]=maxi;
		if(a[j*m+j]==0.0)
      a[j*m+j] = LM_REAL_EPSILON;
		if(j!=m-1){
			tmp=CNST(1.0)/(a[j*m+j]);
			for(i=j+1; i<m; ++i)
        a[i*m+j]*=tmp;
		}
	}

  /* Solve the m linear systems using forward and back substitution */
  for(l=0; l<m; ++l){
      for(i=0; i<m; ++i)  x[i]=0.0;
      x[l]= 1.0;

	    for(i=k=0; i<m; ++i){
		      j=idx[i];
		      sum=x[j];
		      x[j]=x[i];
			    if(k!=0)
			       for(j=k-1; j<i; ++j) sum-=a[i*m+j]*x[j];
		      else
             if(sum!=0.0) k=i+1;
		      x[i]=sum;
	    }

	    for(i=m-1; i>=0; --i){
		      sum=x[i];
		      for(j=i+1; j<m; ++j)
              sum-=a[i*m+j]*x[j];
		          x[i]=sum/a[i*m+i];
	        }

      for(i=0; i<m; ++i) B[i*m+l]=x[i];
  }

  free(buf);
  return 1;
}



/** ==========================================================================
 * Computes in U the uncertainty corresponding to a least squares fit.
 * JtJ:   the approximate Hessian at the solution
 * sumsq: the sum of squared residuals(i.e. goodnes of fit) at the solution,
 * m:     the number of parameters (variables)
 * n:     the number of observations. JtJ can coincide with C.
 * 
 * The function returns the rank of JtJ if successful, 0 on error
 **/
int LM_UNCERTAINTYVAR(LM_REAL *JtJ,
											LM_REAL sumsq,
											int m,
											int n,
											LM_REAL *U) {
   register int i;
   int rnk;
   LM_REAL fact, *C;

   /* temporary matrix */
	 C = (LM_REAL *)malloc( m*m*sizeof(LM_REAL) );
   if(!C){
        fprintf(stderr, "memory allocation request failed\n");
        exit(1);
   }
	
   rnk = LM_LUINVERSE(JtJ, C, m);
   if(!rnk)  return 0;

   rnk = m;

   fact=sumsq/(LM_REAL)(n-rnk);

	 for(i=0; i<m; ++i)
      U[i] = sqrt( C[i*m+i]*fact );
	 
   return rnk;
}



/* ==========================================================================
 * This function seeks the parameter vector p that best describes
 * the measurements vector x. And its uncertainty
 *
 * Returns the number of iterations (>=0) if successfull, -1 if failed
 */

 int LEVMAR_DIF_UN(
  void (*func)(LM_REAL *p, LM_REAL *hx, int m, int n, void *adata),
	/* function describing measurements. A p \in R^m yields a \hat{x} \in  R^n */
  LM_REAL *p,  /* I/O: initial par estimates. Output has estimated solution */
  LM_REAL *x,  /* I: measurement vector */
  int m,       /* I: parameter vector dimension (i.e. #unknowns) */
  int n,       /* I: measurement vector dimension */
  int itmax,   /* I: maximum number of iterations */
  LM_REAL opts[5],
	/* I: opts[0-4] = minim. options [\mu,\epsilon1,\epsilon2,\epsilon3,\delta].
	 */
  LM_REAL info[9],
	/* O: information regarding the minimization. Set to NULL if don't care
	 */
  LM_REAL *work,     /* working memory, allocate if NULL */
  LM_REAL *uncertaintyvar,
	/* O: uncertainty corresponding to LS solution;  Set NULL if not needed. */
  void *adata  /* pointer to possibly additional data */
	)
{
register int i, j, k, l;
int worksz, freework=0, issolved;
/* temp work arrays */
LM_REAL *e,            /* nx1 */
        *hx,           /* \hat{x}_i, nx1 */
        *jacTe,        /* J^T e_i mx1 */
        *jac,          /* nxm */
        *jacTjac,      /* mxm */
        *Dp,           /* mx1 */
        *diag_jacTjac, /* diagonal of J^T J, mx1 */
        *pDp,          /* p + Dp, mx1 */
        *wrk;          /* nx1 */

int using_ffdif=1;
LM_REAL *wrk2=NULL; 
register LM_REAL mu, tmp;
LM_REAL p_eL2,jacTe_inf,pDp_eL2; 
LM_REAL p_L2, Dp_L2=LM_REAL_MAX, dF, dL;
LM_REAL tau, eps1, eps2, eps2_sq, eps3, delta;
LM_REAL init_p_eL2;
int nu, nu2, stop, nfev, njap=0, K=(m>=10)? m: 10, updjac, updp=1, newjac;
const int nm=n*m;

  mu=jacTe_inf=p_L2=0.0;
  stop=updjac=newjac=0;

  if(n<m){
 		 fprintf(stderr, "Can't solve a problem with fewer measurements\n");
     exit(1);
  }

  if(opts){
	  tau     = opts[0];
	  eps1    = opts[1];
	  eps2    = opts[2];
	  eps2_sq = opts[2]*opts[2];
    eps3    = opts[3];
	  delta   = opts[4];
    if(delta<0.0){
       delta=-delta;  /* make positive */
       using_ffdif=0; /* use central differencing */
       wrk2=(LM_REAL *)malloc(n*sizeof(LM_REAL));
       if(!wrk2){
          fprintf(stderr, "Memory allocation request for 'wrk2' failed\n");
          exit(1);
      }
    }
  }
  else{ // use default values
	  tau    = (1.0e-03);
	  eps1   = (1.0e-17);
	  eps2   = (1.0e-17);
	  eps2_sq= (1.0e-34);
    eps3   = (1.0e-17);
	  delta  = (1.0e-06);
  }

  if(!work){
    worksz = LM_DIF_WORKSZ(m, n); //3*n+4*m + n*m + m*m;
    work   =(LM_REAL *)malloc(worksz*sizeof(LM_REAL));
    if(!work){
      fprintf(stderr, "Memory allocation request failed\n");
      exit(1);
    }
    freework=1;
  }

  /* set up work arrays */
  e=work;
  hx=e + n;
  jacTe=hx + n;
  jac=jacTe + m;
  jacTjac=jac + nm;
  Dp=jacTjac + m*m;
  diag_jacTjac=Dp + m;
  pDp=diag_jacTjac + m;
  wrk=pDp + m;

  /* compute e=x - f(p) and its L2 norm */
  (*func)(p, hx, m, n, adata); nfev=1;
  for(i=0, p_eL2=0.0; i<n; ++i){
    e[i]=tmp=x[i]-hx[i];
    p_eL2+=tmp*tmp;
  }
  init_p_eL2=p_eL2;

  nu = 20; /* force computation of J */

  for(k=0; k<itmax; ++k){

    if(p_eL2<=eps3){ /* error is small */
      stop=6;
      break;
    }

    if((updp && nu>16) || updjac==K){/*compute difference approximation to J */
      if(using_ffdif){ /* use forward differences */
        FDIF_FORW_JAC_APPROX(func, p, hx, wrk, delta, jac, m, n, adata);
        ++njap; nfev+=m;
      }
      else{ /* use central differences */
        FDIF_CENT_JAC_APPROX(func, p, wrk, wrk2, delta, jac, m, n, adata);
        ++njap; nfev+=2*m;
      }
      nu=2; updjac=0; updp=0; newjac=1;
    }

    if(newjac){ /* jacobian has changed, recompute J^T J, J^t e, etc */
      newjac=0;

      /* J^T J, J^T e */
      if(nm<=__BLOCKSZ__SQ){      
        for(i=0; i<m; ++i){
          for(j=i; j<m; ++j){
            int lm;

            for(l=0, tmp=0.0; l<n; ++l){
              lm=l*m;
              tmp+=jac[lm+i]*jac[lm+j];
            }

            jacTjac[i*m+j]=jacTjac[j*m+i]=tmp;
          }

          /* J^T e */
          for(l=0, tmp=0.0; l<n; ++l)
            tmp+=jac[l*m+i]*e[l];
          jacTe[i]=tmp;
        }
      }
      else{ 
        /* Cache efficient computation of J^T J based on blocking*/
        TRANS_MAT_MAT_MULT(jac, jacTjac, n, m, __BLOCKSZ__);

        /* cache efficient computation of J^T e */
        for(i=0; i<m; ++i)
          jacTe[i]=0.0;

        for(i=0; i<n; ++i){
          register LM_REAL *jacrow;

          for(l=0, jacrow=jac+i*m, tmp=e[i]; l<m; ++l)
            jacTe[l]+=jacrow[l]*tmp;
        }
      }
      
      /* Compute ||J^T e||_inf and ||p||^2 */
      for(i=0, p_L2=jacTe_inf=0.0; i<m; ++i){
         if(jacTe_inf < (tmp=FABS(jacTe[i]))) jacTe_inf=tmp;

         diag_jacTjac[i]=jacTjac[i*m+i]; 
         p_L2+=p[i]*p[i];
      }
    }

    /* check for convergence */
    if((jacTe_inf <= eps1)){
      Dp_L2=0.0; /* no increment for p in this case */
      stop=1;
      break;
    }

   /* compute initial damping factor */
    if(k==0){
      for(i=0, tmp=LM_REAL_MIN; i<m; ++i)
        if(diag_jacTjac[i]>tmp)
					 tmp=diag_jacTjac[i]; /* find max diagonal element */
      mu=tau*tmp;
    }

    /* determine increment using adaptive damping */

    /* augment normal equations */
    for(i=0; i<m; ++i)  jacTjac[i*m+i]+=mu;

    /* solve augmented equations */
    /* use the LU included with levmar */
    issolved = AX_EQ_B_LU(jacTjac, jacTe, Dp, m);

    if(issolved){
    /* compute p's new estimate and ||Dp||^2 */
      for(i=0, Dp_L2=0.0; i<m; ++i){
          pDp[i]=p[i] + (tmp=Dp[i]);
          Dp_L2+=tmp*tmp;
      }

      if(Dp_L2<=eps2_sq*p_L2){ /* relative change in p is small, stop */
        stop=2;
        break;
      }

      if(Dp_L2>=(p_L2+eps2)/(CNST(EPSILON)*CNST(EPSILON))){/*almost singular */
        stop=4;
        break;
      }

      (*func)(pDp, wrk, m, n, adata); ++nfev; /* evaluate function at p + Dp */
      for(i=0, pDp_eL2=0.0; i<n; ++i){ /* compute ||e(pDp)||_2 */
        tmp=x[i]-wrk[i];
        pDp_eL2+=tmp*tmp;
      }

      dF=p_eL2-pDp_eL2;
      if(updp || dF>0){ /* update jac */
        for(i=0; i<n; ++i){
          for(l=0, tmp=0.0; l<m; ++l)
            tmp+=jac[i*m+l]*Dp[l]; /* (J * Dp)[i] */
          tmp=(wrk[i] - hx[i] - tmp)/Dp_L2;
          for(j=0; j<m; ++j)
            jac[i*m+j]+=tmp*Dp[j];
        }
        ++updjac;
        newjac=1;
      }

      for(i=0, dL=0.0; i<m; ++i)
         dL+=Dp[i]*(mu*Dp[i]+jacTe[i]);

      if(dL>0.0 && dF>0.0){ /* reduction in error, increment is accepted */
         dF=(CNST(2.0)*dF/dL-CNST(1.0));
         tmp=dF*dF*dF;
         tmp=CNST(1.0)-tmp*tmp*dF;
         mu=mu*( (tmp>=CNST(ONE_THIRD))? tmp : CNST(ONE_THIRD) );
         nu=2;

         for(i=0 ; i<m; ++i) /* update p's estimate */
           p[i]=pDp[i];

         for(i=0; i<n; ++i){ /* update e, hx and ||e||_2 */
           e[i]=x[i]-wrk[i];
           hx[i]=wrk[i];
         }
         p_eL2=pDp_eL2;
         updp=1;
         continue;
      }
    }


    mu *= nu;
    nu2 = nu<<1; 
    if(nu2<=nu){ /* nu has wrapped around (overflown). */
      stop=5;
      break;
    }
    nu=nu2;

    for(i=0; i<m; ++i) /* restore diagonal J^T J entries */
      jacTjac[i*m+i]=diag_jacTjac[i];
  }

  if(k>=itmax) stop=3;

  for(i=0; i<m; ++i) /* restore diagonal J^T J entries */
    jacTjac[i*m+i]=diag_jacTjac[i];

  if(info){
    info[0] = init_p_eL2;
    info[1] = p_eL2;
    info[2] = jacTe_inf;
    info[3] = Dp_L2;
    for(i=0, tmp=LM_REAL_MIN; i<m; ++i)
       if(tmp<jacTjac[i*m+i]) tmp=jacTjac[i*m+i];
    info[4] = mu/tmp;
    info[5] = (LM_REAL)k;
    info[6] = (LM_REAL)stop;
    info[7] = (LM_REAL)nfev;
    info[8] = (LM_REAL)njap;
  }

  /* uncertainty */
  if(uncertaintyvar){
     LM_UNCERTAINTYVAR(jacTjac, p_eL2, m, n, uncertaintyvar);
  }

                                                               
  if(freework) free(work);

  if(wrk2) free(wrk2);

  return (stop!=4)?  k : -1;
}

/* undefine everything. THIS MUST REMAIN AT THE END OF THE FILE */
#undef LEVMAR_DIF_UN
#undef FDIF_FORW_JAC_APPROX
#undef FDIF_CENT_JAC_APPROX
#undef TRANS_MAT_MAT_MULT
#undef AX_EQ_B_LU
