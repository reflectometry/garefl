/* This program is public domain. */

#ifndef _INTERFACE_H
#define _INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

/* ==== interface handling ==== */

/* Interface handling:
 *
 * To create a new interface style, define a function
 *   f(n,dz,r)
 * which fills z with the width of each interface slice
 * and fills r with a value from -1 to 1, where 0 means
 * average the values on either side of the interface,
 * -1 means favour the top of the interface and 1 means
 * favour the bottom of the interface.  The ideal is to
 * choose z=finv(r) so that the steps in r are uniform from
 * -1+1/2n by 1/n to 1-1/2n, then convert this to slices
 * using interface_steps(n,z).
 *
 * If n == 0, set *r to f(*dz) --- we need this for the
 * case where we are not choosing the steps in the interface
 * but instead are forced (due to overlapping interfaces
 * or due to functional layers) to have interfaces at
 * particular positions.  [We may later decide to interpolate
 * rather than use this functionality].
 *
 * tanh_interface
 * erf_interface
 *   Predefined interface functions
 * interface_steps(int n, Real z[])
 *   convert from depths at z to widths of slices whose edges fall between z
 * interface_create(interface*, const char *name, interface_fn fn, int n)
 *   create an interface vector and populate it with n steps.  Name
 *   must refer to persistent memory.
 * interface_destroy(interface*)
 *   free memory for the interface vector
 * Real interface_value(interface*, Real z)
 *   compute the value of the interface at z (either by interpolation
 *   or by calling to the interface function directly if that is not
 *   too slow).
 * Real interface_average(Real v1, v2, r)
 *   given an interface value f(z)=r and values v1 and v2 on either
 *   side of the interface, compute the diffuse value at z.
 * Real interface_overlap(Real r1, r2, v1, v, v2)
 *   given overlapping interfaces at point z, we have z1 as the distance
 *   from the first interface to z and z2 as the distance from the second
 *   interface to z.  With f1(z1)=r1 and f2(z2)=r2, and with v1 the value
 *   before the first interface, v the nominal value between the interfaces
 *   and v2 the value at the second interface, compute the diffuse value at z.
 */
typedef void (*interface_fn)(int n, Real dz[], Real r[]);
typedef struct interface_struct {
  int n;
  Real depth;
  Real *value;
  Real *step;
  interface_fn fn;
  const char *name;
} interface;
void tanh_interface(int n, Real dz[], Real r[]);
void erf_interface(int n, Real dz[], Real r[]);
void interface_steps(int n, Real dz[]);
void interface_init(interface *rm);
int interface_create(interface *rm, const char *name,
		     interface_fn f, int n);
int interface_set(interface *rm, const char *name,
		  int n, const Real z[], const Real v[]);
void interface_destroy(interface *rm);
Real interface_value(const interface *rm, Real v);
Real interface_average(Real above, Real below, Real weight);
Real interface_overlap(Real weight_above, Real weight_below,
			 Real v_above, Real v, Real v_below);
void interface_print(interface *fm);

#ifdef __cplusplus
}
#endif

#endif

/* $Id$ */
