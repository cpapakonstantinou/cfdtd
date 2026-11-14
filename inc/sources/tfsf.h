/**
 * @file sources/tfsf.h
 * @brief TF/SF corrections of the plane-wave source on the six faces of the TF/SF box
 *
 * The incident field comes from the 1D auxiliary grid, sources/aux_grid.h. On each face the H (or E)
 * cells just outside, inside, the box get the incident E (H) of the 1D grid times the update
 * coefficient, so the total field lives inside the box and the scattered field outside.
 * @verbatim
      face  H corrected                         E corrected
      x     Hy, from Ez_inc, pol z,             Ey, from Hz_inc, pol x
      y     Hx (Ez_inc), Hz (Ex_inc)            Ex (Hz_inc), Ez (Hx_inc)
      z     Hy, Ex_inc, pol x,                  Ey, Hx_inc, pol z
   @endverbatim
 * The near y face also carries the mode source, mode profile times the signal, and, with periodic x
 * and z, takes the incident field from the reference line.
 *
 * E coefficient of the background m: ce_2, dielectric, ce_2_d, Drude, Lorentz; the other dispersive
 * models add the TF/SF terms inside their own update, tfsf_e_terms, handed to the boundary as its incident field: boundary/incident_field.h.
 */
#ifndef FDTD_SOURCES_TFSF_H
#define FDTD_SOURCES_TFSF_H

namespace fdtd::sources
{
enum class Face { XNea, XFar, YNea, YFar, ZNea, ZFar };

/// H correction on one face of this rank
void tfsf_h(Face f);
/// E correction on one face of this rank
void tfsf_e(Face f);

/**
 * @brief The TF/SF terms of E component A at cell (i, j, k), scaled by c, for the dispersive models
 *        that add them inside their update, they are not in the background rule of tfsf_e
 * @verbatim
      Ex: y planes, pol x,      Ey: z planes, pol z, then x planes, pol x,      Ez: y planes, pol z
   @endverbatim
 */
void tfsf_e_terms(int A, int i, int j, int k, double c, double &e_new);
} // namespace fdtd::sources

#endif
