/** \file   core/limits.h
 *  \brief  Capacity limits of the solver's fixed-size tables.
 */
#ifndef FDTD_CORE_LIMITS_H
#define FDTD_CORE_LIMITS_H

namespace fdtd
{
inline constexpr int max_frequencies = 401; ///< analysed, DFT, frequencies, at most: size of the per-frequency spectra tables
} // namespace fdtd

#endif
