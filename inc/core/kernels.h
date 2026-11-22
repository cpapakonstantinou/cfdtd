/**
 * @file core/kernels.h
 * @brief Contiguous vector kernels used by the update loops.
 *
 * Every kernel works on contiguous runs, pointer + length. The default backend is a plain loop
 * that the compiler vectorises at -O2/-O3. Building with `make DAMM=<path to damm>` defines
 * FDTD_DAMM and routes the kernels through DAMM, https://github.com/cpapakonstantinou/damm;
 * that needs DAMM's `fused_union` remainder fix, docs/dev/damm-fused_union-rem_cols.patch.
 *
 * The plain backend does the same arithmetic per element as the scalar loops it replaces, one
 * multiply, one add, so `make MARCH=` builds are bit-identical to them. DAMM fuses the pair into an
 * FMA and needs -march=native: results differ by round-off, at most 9e-12 relative on the regression
 * case, volume DFT and exact multipoles included.
 */
#ifndef FDTD_KERNELS_H
#define FDTD_KERNELS_H

#include <complex>
#include <cstddef>

#ifdef FDTD_DAMM
#include <damm.h>
#endif

namespace kernels
{

/**
 * @brief Running-DFT step over n contiguous values: re += x u, im += x v
 * @verbatim
     x  : | x0 | x1 | x2 | ... | x(n-1) |     one plane of a field, gathered once
     re : | r0 | r1 | r2 | ... |        |  += x * u    u = cos(w t)
     im : | i0 | i1 | i2 | ... |        |  += x * v    v = sin(w t)
   @endverbatim
 */
template <typename T>
inline void dft_accumulate(T *__restrict re, T *__restrict im, T u, T v, const T *__restrict x, std::size_t n)
{
#ifdef FDTD_DAMM
	T *R[1] = {re}, *I[1] = {im}, *X[1] = {const_cast<T *>(x)};
	using damm::FusionPolicy;
	damm::scalar::fused_union<FusionPolicy::FUSION_FIRST, T, std::plus<>, std::multiplies<>>(R, u, X, R, 1, n);
	damm::scalar::fused_union<FusionPolicy::FUSION_FIRST, T, std::plus<>, std::multiplies<>>(I, v, X, I, 1, n);
#else
	for (std::size_t q = 0; q < n; ++q)
	{
		re[q] += x[q] * u;
		im[q] += x[q] * v;
	}
#endif
}

/// c[q] += x[q] (u + i v): the real and imaginary parts as dft_accumulate, on one complex row
template <typename T>
inline void dft_accumulate(std::complex<T> *__restrict c, T u, T v, const T *__restrict x, std::size_t n)
{
	T *p = reinterpret_cast<T *>(c);
	for (std::size_t q = 0; q < n; ++q)
	{
		p[2 * q] += x[q] * u;
		p[2 * q + 1] += x[q] * v;
	}
}

} // namespace kernels

#endif // FDTD_KERNELS_H
