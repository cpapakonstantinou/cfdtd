/**
 * @file materials/models/cell.h
 * @brief What a material model's point update receives, and the shape every model has.
 *
 * A model is a type with one static member template
 * @code
 *   template <int A, class Boundary>
 *   static void update(const Cell &c, Boundary &&boundary);
 * @endcode
 * that advances E component A (0 x, 1 y, 2 z) and the model's own state at one cell. Its arrays come
 * from materials::Component<A>; the parts that differ per component are written with the cyclic axes
 * B1 = A+1, B2 = A+2. `boundary(coeff, e_new)` adds the CPML psi and TF/SF terms scaled by the model's
 * curl coefficient. Drude and Lorentz do not call it: their CPML psi and TF/SF terms are added after
 * the scan, boundary::e_face, planeWave.cc, scaled by ce_2_d.
 *
 * The scan, src/update/update_e.cc, picks the model once per stretch, materials::with_model, and calls
 * update<A> for its cells.
 */
#ifndef FDTD_MATERIALS_MODELS_CELL_H
#define FDTD_MATERIALS_MODELS_CELL_H

namespace fdtd::materials::models
{
/// one cell of a stretch
struct Cell
{
	int mm;         ///< dispersive material index, code - DispMateOffsetD
	int i, j, k;    ///< local cell
	int J2;         ///< global y index, nonlocal barrier
	double curl;    ///< curl H of the field, CPML-stretched outside the inner region
	double curl_nl; ///< curl H of the harmonic field, separate linear / nonlinear fields, else 0
};
} // namespace fdtd::materials::models

#endif
