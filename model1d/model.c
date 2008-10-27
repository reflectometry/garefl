/* This program is public domain. */

/* fzw */
/** \file
 *  Handle model.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "model.h"

/* The usual model size, if no size is given; we dynamically resize as
 * needed, so it is not so important that this number be accurate. */
#define MODEL_LENGTH_DEFAULT 30


/* Subroutine generates layer profile for the refractive indices QCSQ,
   absorptions MU, and interfaces ROUGH specified
   by the fit program by explicitly calculating the interfaces of the
   interfacial regions as variations of average QCSQ and MU on a hyperbolic
   tangent or error function profile [see Anastasiadis, Russell, Satija, and
   Majkrzak, J. Chem. Phys. 92, 5677 (1990)]
   John Ankner 14-June-1990 */

/* FIXME this algorithm is a bit silly: it uses an equal number of interface 
 * steps for each interface no matter how wide the interface or how
 * big the difference in refractive index.
 */

static const char* errmsg[] = {
  /*   0 */ "unknown model error",
  /*  -1 */ "memory allocation error",
  /*  -2 */ "too many repeats for model",
  /*  -3 */ "bad repeat parameters",
  /*  -4 */ "repeat overlaps with other repeats",
} ;

const char *model_error(int code)
{
  code = -code;
  if (code < 0 || code >= sizeof(errmsg)/sizeof(errmsg[0])) code = 0;
  return errmsg[code];
}

void model_print(const model *m, const char *filename)
{
  int i;
  FILE *f;

  if (filename == NULL) {
    f = stdout;
  } else {
    f = fopen(filename, "w");
    if (f == NULL) return;
  }
    
  fprintf(f,"Model %p\n",m);
  fprintf(f,"  size %d out of %d\n",m->n,m->capacity);
  fprintf(f,"  interface=%p\n",m->rm);
  fprintf(f,"  repeats=%d at %p\n",m->num_repeats,m->repeats);
  for (i=0; i< m->num_repeats; i++) {
    fprintf(f,"    %d %d %d %d\n",i,
	    m->repeats[3*i],m->repeats[3*i+1],m->repeats[3*i+2]);
  }
#ifdef HAVE_MAGNETIC
  fprintf(f,"  # d rho mu rough P Prough theta thetarough\n");
  for (i=0; i < m->n; i++) {
    fprintf(f,"  %d %g %g %g %g %g %g %g %g\n", i, m->d[i], 
	   m->rho[i], m->mu[i], m->rough[i], 
	   m->P[i], m->Prough[i], m->theta[i], m->thetarough[i]);
  }
#else
  fprintf(f,"  # d(%p) rho(%p) mu(%p) rough(%p)\n",m->d,m->rho,m->mu,m->rough);
  for (i=0; i < m->n; i++) {
    fprintf(f,"  %d %g %g %g %g\n",i,m->d[i],m->rho[i],m->mu[i],m->rough[i]);
  }
#endif

  if (f != stdout) fclose(f);
}


void model_init(model *m)
{
  m->num_repeats=0;
  m->capacity = 0;
  m->n = 0;
  m->is_magnetic = 0;
}

int model_extend(model *m, int n)
{
  int i;

  if (m->capacity < 0) return 0;
  if (m->n+n > m->capacity) {
    int old_size = m->capacity;
    int new_size = m->capacity + n;
    new_size += new_size/10 + 20; /* 10% spare */
    if (new_size < MODEL_LENGTH_DEFAULT) new_size = MODEL_LENGTH_DEFAULT;
    if (m->capacity <= 0) {
      m->mu = malloc(MODEL_FIELDS*sizeof(double)*new_size);
    } else {
      m->mu = realloc(m->mu,MODEL_FIELDS*sizeof(double)*new_size);
      if (m->mu != NULL) {
	for (i=MODEL_FIELDS-1; i>0; i--) {
	  memmove(m->mu+i*new_size,m->mu+i*old_size,old_size*sizeof(double));
	}
      }
    }
    m->rho = m->mu + new_size;
    m->d = m->mu + 2*new_size;
    m->rough = m->mu + 3*new_size;
#ifdef HAVE_MAGNETIC
    m->P = m->mu + 4*new_size;
    m->Prough = m->mu + 5*new_size;
    m->theta = m->mu + 6*new_size;
    m->thetarough = m->mu + 7*new_size;
#endif
    m->capacity = new_size;
    if (m->mu==NULL) {
      model_destroy(m);
      m->capacity = -1; /* remember that we destroyed the profile */
      return 0;
    }
  }
  return 1;
}

void model_destroy(model *m)
{
  if (m->capacity > 0) {
    if (m->mu != NULL) free(m->mu);
  }
  model_init(m);
}

void 
model_layer(model *m, double d, double rho, double mu, double rough)
{
  int n = m->n;

  model_extend(m,1);
  if (m->n >= m->capacity) return;
#ifdef HAVE_MAGNETIC
  m->P[n]=0;
  m->Prough[n]=rough;
  m->theta[n]=270.;
  m->thetarough[n]=rough;
#endif
  m->rho[n]=rho;
  m->mu[n]=mu;
  m->d[n] = d;
  m->rough[n]=rough;
  m->n=n+1;
}

void 
model_magnetic(model *m, double d, double rho, double mu, double rough,
	       double P, double Prough, double theta, double thetarough)
{
#ifdef HAVE_MAGNETIC
  int n = m->n;

  m->is_magnetic = 1;
  model_extend(m,1);
  if (m->n >= m->capacity) return;
  m->P[n]=P;
  m->Prough[n]=Prough;
  m->theta[n]=theta;
  m->thetarough[n]=thetarough;
  m->rho[n]=rho;
  m->mu[n]=mu;
  m->d[n] = d;
  m->rough[n]=rough;
  m->n=n+1;
#else
  fprintf(stderr,"librefl not compiled for magnetic systems\n");
  exit(1);
#endif
}


/* ================================================================== */
/* Repeat handling */

/* Assumptions:
 * repeats are ordered
 * repeats don't overlap
 * there is a maximum number of repeats
 * repeats is a fittable parameter, constrained to be an integer
 * repeat_pars in model must point somewhere in the parameter list
 * there are no repeats of length 1
 */
void model_repeat(model *m, int R, int *R_start, int *R_end, int *R_count)
{
  if (R > m->num_repeats) {
    /* no more repeats --- specify a set of layers outside the model */
    *R_start = -1;
    *R_end = -1;
    *R_count = 0;
  } else {
    *R_start = m->repeats[3*R];
    *R_end = m->repeats[3*R+1];
    *R_count = m->repeats[3*R+3];
  }
}

int model_repeat_insert(model *m, int R_start, int R_end, int R_count)
{
  int r, pos;
  if (m->num_repeats >= MODEL_MAX_REPEATS) return MODEL_REPEATS_EXCEEDED;
  if (R_end < R_start+1 || R_count < 2) return MODEL_REPEATS_BAD;
  for (pos=0; pos < m->num_repeats; pos++) 
    if (m->repeats[2*pos] < R_start) break;
  if ( (pos > 0 && m->repeats[2*pos-1]>=R_start)
       || (pos < m->num_repeats && m->repeats[2*pos]>=R_end) )
    return MODEL_REPEATS_OVERLAP;

  /* Shift tail of repeat structure out of the way */
  m->num_repeats++;
  for (r=m->num_repeats; r>pos; r--) {
    m->repeats[3*r]=m->repeats[3*r-3];
    m->repeats[3*r+1]=m->repeats[3*r-2];
    m->repeats[3*r+2]=m->repeats[3*r-1];
  }

  m->num_repeats++;
  m->repeats[3*pos] = R_start;
  m->repeats[3*pos+1] = R_end;
  m->repeats[3*pos+2] = R_count;

  return 0;
}

void model_repeat_delete(model *m, int R)
{
  int r;
  for (r=R+1; r < m->num_repeats; r++) {
    m->repeats[3*r-3]=m->repeats[3*r];
    m->repeats[3*r-2]=m->repeats[3*r+1];
    m->repeats[3*r-1]=m->repeats[3*r+2];
  }
  m->num_repeats--;
}



/* ================================================================= */
/* Profile generation */


/* Note: midpoint is chosen so that the central interface slice is
 * placed at the end of the layer rather than at the beginning when
 * the number of interface slices is odd. */

#ifdef HAVE_MAGNETIC
/* This is the minimum width at which we combine roughness */
#define MINSTEP 0.01
static void
add_half_interface(model *m, profile *p, int above, int below, int leftside)
{
  interface *rm=m->rm;
  int i, j, midpoint = (rm->n+1)/2;
  double rP = m->Prough[below];
  double rrho = m->rough[below];
  double rtheta = m->thetarough[below];
  double sm, med, big, z, c;
  int n = p->n;


#ifdef VERBOSE
  printf("Add %s interface of width %g:\n\trho %g -> %g\n",
	 leftside ? "top" : "bottom",
	 m->rough[below], m->rho[above], m->rho[below]);
#endif

  /* At worst, we will have 3 times the number of slices in one interface */
  if (!profile_extend(p,3*midpoint)) return;

  /* Sort interface widths as sm, med, big */
  if (rP < rrho)
    if (rP > rtheta) sm=rtheta, med=rP, big = rrho;
    else if (rrho < rtheta) sm=rP, med=rrho, big=rtheta;
    else sm=rP, med=rtheta, big=rrho;
  else /* rP > rrho */
    if (rP < rtheta) sm=rrho, med=rP, big=rtheta;
    else if (rrho < rtheta) sm=rrho, med=rtheta, big=rP;
    else sm=rtheta, med=rrho, big=rP;

  /* If the interface has no width, do nothing */
  if (big < 1e-10) return;

  /* We are going to first add the slices for the sharpest interface,
   * then extend it with slices remaining from the next sharpest
   * interface, and finally extend it with slices from the broadest
   * interface. This will make sure that we are sampling as densely
   * as needed for the steepest curve. 
   *
   * We are storing the widths of the slices in the profile slabs that
   * we have reserved above but not yet used.  For the bottom side
   * interface, we store the slices upsidedown, from the center of
   * the interface to the edge, then reverse them.
   */
  z = 0.;
  if (leftside) {
    // printf("left side %d->%d\n",above,below);

    /* add sharpest interface (reverse order) */
    for (i=midpoint-1; i >=0; i--) z += (p->d[n++] = sm*rm->step[i]);
    // printf("sm: i=%d out of %d, d=%g, z=%g, n=%d\n",i,rm->n,p->d[n-1],z,n);

    /* add extra slices to fill out the medium interface (reverse order) */
    c = 0.;
    for (i=midpoint-1; i>=0 && c < z; i--) c += med*rm->step[i];
    // printf("med overlap i=%d, c=%g, z=%G ?\n",i,c,z);
    if (c-z < MINSTEP) p->d[n-1] += c-z;
    else p->d[n++] = c-z;
    z = c;
    while (i >= 0) z += (p->d[n++] = med*rm->step[i--]);
    // printf("med extend: i=%d out of %d, z=%g, n=%d\n",i,rm->n,z,n);

    /* add extra slices to fill out the widest interface (reverse order) */
    c = 0.;
    for (i=midpoint-1; i>=0 && c < z; i--) c += big*rm->step[i];
    // printf("big overlap i=%d, c=%g, z=%G ?\n",i,c,z);
    if (c-z < MINSTEP) p->d[n-1] += c-z;
    else p->d[n++] = c-z;
    z = c;
    while (i >= 0) z += (p->d[n++] = big*rm->step[i--]);
    // printf("big extend: i=%d out of %d, z=%g, n=%d\n",i,rm->n,z,n);

    /* reverse the slices, and measure from -z to 0 rather than 0 to z */
    // printf("reversing\n");
    for (i=p->n, j=n-1; i < j; i++,j--) {
      double temp = p->d[i];
      p->d[i] = p->d[j];
      p->d[j] = temp;
    }
    z = -z;
  }
  else {
    // printf("right side %d->%d\n",above,below);

    /* add sharpest interface */
    for (i=midpoint; i < rm->n; i++) z += (p->d[n++] = sm*rm->step[i]);
    // printf("sm: i=%d out of %d, d=%g, z=%g, n=%d, sm=%g, step=%g\n",i,rm->n,p->d[n-1],z,n,sm,rm->step[i]);

    /* add extra slices to fill out the medium interface */
    c = 0.;
    for (i=midpoint; i < rm->n && c < z; i++) c += med*rm->step[i];
    // printf("med overlap i=%d out of %d, c=%g, n=%d\n",i,rm->n,c,n);
    if (c-z < MINSTEP) p->d[n-1] += c-z;
    else p->d[n++] = c-z;
    z = c;
    while (i < rm->n) z += (p->d[n++] = med*rm->step[i++]);
    // printf("med extend: i=%d out of %d, z=%g, n=%d\n",i,rm->n,z,n);

    /* add remaining slices to fill out the widest interface */
    c = 0.;
    for (i=midpoint; i < rm->n && c < z; i++) c += big*rm->step[i];
    // printf("big overlap i=%d out of %d, c=%g, n=%d\n",i,rm->n,c,n);
    if (c-z < MINSTEP) p->d[n-1] += c-z;
    else p->d[n++] = c-z;
    z = c;
    while (i < rm->n) z += (p->d[n++] = big*rm->step[i++]);
    // printf("big extend: i=%d out of %d, c=%g, n=%d\n",i,rm->n,c,n);

    z = 0.;
  }

  /* Fill in the slices */
  while (p->n < n) {
    double d = p->d[p->n];
    double zmid = z + d/2;
    double wrho = interface_value(m->rm, zmid/rrho);
    double wtheta = interface_value(m->rm, zmid/rtheta);
    double wP = interface_value(m->rm, zmid/rP);
    // printf("add_half_interface slice %d: Adding %g at %g using %d->%d with %g\n",p->n,d,zmid,above,below,wrho);
    profile_slice(p, d,
		  interface_average(m->rho[above],m->rho[below],wrho),
		  interface_average(m->mu[above],m->mu[below],wrho),
		  interface_average(m->P[above],m->P[below],wP),
		  interface_average(m->theta[above],m->theta[below],wtheta));
    z += d;
  }
}

static void
add_interface_left(model *m, profile *p, int above, int below)
{
  add_half_interface(m,p,above,below,1);
}

static void
add_interface_right(model *m, profile *p, int above, int below)
{
  add_half_interface(m,p,above,below,0);
}

#else

static void
add_interface_left(model *m, profile *p, int above, int below)
{
  interface *rm=m->rm;
  int i, midpoint = (rm->n+1)/2;
 
  if (!profile_extend(p,midpoint)) return;

  /* If the interface has no width, do nothing */
  if (m->rough[below] < 1e-10) return;

  for (i = 0; i < midpoint; i++) {
    // printf("in add_interface_left\n");
    profile_slice(p, rm->step[i]*m->rough[below],
		  interface_average(m->rho[above],m->rho[below],rm->value[i]),
		  interface_average(m->mu[above],m->mu[below],rm->value[i]));
  }
}

static void
add_interface_right(model *m, profile *p, int above, int below)
{
  interface *rm=m->rm;
  int i, midpoint = (rm->n+1)/2;
 
  if (!profile_extend(p,rm->n-midpoint+1)) return;

  /* If the interface has no width, do nothing */
  if (m->rough[above] < 1e-10) return;

  for (i = midpoint; i < rm->n; i++) {
    // printf("in add_interface_right\n");
    profile_slice(p, rm->step[i]*m->rough[below],
		  interface_average(m->rho[above],m->rho[below],rm->value[i]),
		  interface_average(m->mu[above],m->mu[below],rm->value[i]));
  }
}
#endif


static void
add_overlapping_slice(model *m, profile *p, int above, int layer, int below,
		      double position, double width)
{
  double rho_weight_above =
    interface_value(m->rm, position/m->rough[layer]);
  double rho_weight_below =
    interface_value(m->rm, (position-m->d[layer])/m->rough[below]);
#ifdef HAVE_MAGNETIC
  double P_weight_above =
    interface_value(m->rm, position/m->Prough[layer]);
  double P_weight_below =
    interface_value(m->rm, (position-m->d[layer])/m->Prough[below]);
  double theta_weight_above =
    interface_value(m->rm, position/m->thetarough[layer]);
  double theta_weight_below =
    interface_value(m->rm, (position-m->d[layer])/m->thetarough[below]);
#endif

  // printf("in add_overlapping_slice\n");
  profile_slice(p, width,
		interface_overlap(rho_weight_above, rho_weight_below,
				  m->rho[above],m->rho[layer],m->rho[below]),
		interface_overlap(rho_weight_above, rho_weight_below,
				  m->mu[above],m->mu[layer],m->mu[below])
#ifdef HAVE_MAGNETIC
		,interface_overlap(P_weight_above, P_weight_below,
				  m->P[above],m->P[layer],m->P[below])
		,interface_overlap(theta_weight_above, theta_weight_below,
			  m->theta[above],m->theta[layer],m->theta[below])
#endif
		);
}

static void 
add_layer(model *m, profile *p, int above, int layer, int below)
{
  double rtop, rbottom, bulk_depth;
  interface *rm = m->rm;
  /* Make sure there is room for the slices, and then some */
  if (!profile_extend(p,2*rm->n+2)) return;

  /* Find non-overlapping portion */
  rtop = m->rough[layer];
  rbottom = m->rough[below];
#ifdef HAVE_MAGNETIC
  if (m->Prough[layer] > rtop) rtop=m->Prough[layer];
  if (m->Prough[below] > rbottom) rbottom=m->Prough[below];
  if (m->thetarough[layer] > rtop) rtop=m->thetarough[layer];
  if (m->thetarough[below] > rbottom) rbottom=m->thetarough[below];
#endif
  bulk_depth = m->d[layer] - rm->depth*(rtop+rbottom)/2.;
  if (bulk_depth > 1e-10) {
    /* Interfaces don't overlap, so add top and bottom interfaces
     * separately, with a bulk layer in between. */
    add_interface_right(m,p,above,layer);
    // printf("in add_layer with bulk_depth=%g\n",bulk_depth);
    profile_slice(p,bulk_depth,m->rho[layer],m->mu[layer]
#ifdef HAVE_MAGNETIC
		  , m->P[layer], m->theta[layer]
#endif
		  );
    add_interface_left(m,p,layer,below);
  } else {
    /* Interfaces overlap */
    int n = rm->n;
#ifdef VERBOSE
    printf("Add overlapping interface of width %g:\n\trho %g -> %g -> %g\n\tmu %g -> %g -> %g\n",
	   m->d[above],
	   m->rho[above], m->rho[layer], m->rho[below],
	   m->mu[above], m->mu[layer], m->mu[below]);
#endif
#if 0
    /* Simple model: equal width slices sampled at the centers */
    double step = m->d[layer] / (n+1.);
    int i;

    for (i=0; i <= n; i++)
      add_overlapping_slice(m,p,above,layer,below,step*(i+.5),step);
#else
    /* Follow the algorithm in the mlayer program, where the first 
     * and last slices are half-width and the remaining equally-spaced 
     * slices are full-width.
     */
    double step = m->d[layer] / n;
    int i;

    add_overlapping_slice(m,p,above,layer,below,step/4.,step/2.);
    for (i=1; i < n; i++)
      add_overlapping_slice(m,p,above,layer,below,step*i,step);
    add_overlapping_slice(m,p,above,layer,below,m->d[layer]-step/4.,step/2.);
#endif
  }
}


void model_profile(model *m, profile *p)
{
  int layer, repeat_start, R, R_start, R_end, R_count;

  profile_reset(p);

  /* Insert vacuum layers */
  if (!profile_extend(p,1)) return;
  profile_slice(p,10.,m->rho[0],m->mu[0]
#ifdef HAVE_MAGNETIC
		,m->P[0],m->theta[0]
#endif
		);
  add_interface_left(m,p,0,1);
  p->vacuum_offset = profile_depth(p); /* vacuum interface thickness */

  /* initialize R_start, R_end, R_count */
  model_repeat(m,R=0,&R_start,&R_end,&R_count);

  repeat_start = -1;
  for (layer=1; layer<m->n-1; layer++) {
    /* Remember where the repeat section starts */
    if (layer == R_start+1) repeat_start = p->n;
    
    /* Check if we need to repeat */
    if (layer == R_end) {
      /* The interface between repeats is different from the
       * interface into the first repeat and out of the last
       * repeat. So instead of repeating R1 R2 R3 R4 as in:
       *       Ln|R1 R2 R3 R4|R1 R2 R3 R4|R1 R2 R3 R4|Ln+1
       * we have to repeat R2 R3 R1 as in:
       *       Ln R1*|R2 R3 R4 R1|R2 R3 R4 R1|R2 R3 R4* Ln+1
       * with the last repeat only extending from R_start+1 to
       * R_end-1.
       *
       * At this point we have:
       *       Ln R1*|R2 R3
       * so we need to add R4 and R1 to complete the repeating
       * section, then copy that section R_count-2 times (first
       * repeat is done, last repeat is special), then copy
       * what we have of the last repeat, leaving us at:
       *       Ln R1* R2 R3 R4 R1 R2 R3 R4 R1 R2 R3
       * and layer == R4, so continue as normal. */
      int r, repeat_length, repeat_interface;

      repeat_interface = p->n;
      
      /* add interface between repeating sections */
      add_layer(m,p,R_end-1,R_end,R_start);
      add_layer(m,p,R_end,R_start,R_start+1);
      
      /* do bulk repeats of R_start+1 to R_start */
      assert(repeat_start > 0);
      repeat_length = p->n - repeat_start;
      profile_extend(p,R_count*repeat_length); /* make sure there is room */
      for (r=1; r < R_count-1; r++)
	profile_copy(p, repeat_start, repeat_length);

      /* copy layers from R_start to I_start */
      profile_copy(p, repeat_start, repeat_interface-repeat_start);

      /* Look for next repeat */
      model_repeat(m,++R,&R_start,&R_end,&R_count);
      repeat_start = -1;
    }
    add_layer(m,p,layer-1,layer,layer+1);
  }

  /* Insert substrate layers */
  add_interface_right(m,p,layer-1,layer);
  if (!profile_extend(p,1)) return;
  profile_slice(p,10.,m->rho[layer],m->mu[layer]
#ifdef HAVE_MAGNETIC
		,m->P[layer],m->theta[layer]
#endif
		);

  /* Convert polar to cartesian for magnetic layers */
  profile_expth(p);
}

/* $Id$ */
