/**
 * @file sources/mode_spec.h
 * @brief The mode source's parameters, [mode] in the configuration
 */
#ifndef FDTD_SOURCES_MODE_SPEC_H
#define FDTD_SOURCES_MODE_SPEC_H

#include <string>

namespace fdtd::sources
{
/** [mode]: a waveguide mode from MPB, injected one way on a plane, [source] type = mode.
 *  The profile is read from MPB's field files of one band, e.k01.b01.h5, h.k01.b01.h5: a 2D
 *  cross-section, computed with k along the waveguide. MPB's cell is centred on the structure
 *  centre, plus `offset_*`, its x and y axes lie along `mpb_x_axis` and `mpb_y_axis`. */
struct ModeSpec
{
	std::string efield;                 ///< MPB E-field file, relative to [paths] input_dir
	std::string hfield;                 ///< MPB H-field file
	double unit_length = 1e-6;          ///< [m] MPB's length unit a
	double effective_index = 0;         ///< n_eff = k / frequency of the band, MPB units
	std::string propagation_axis = "y"; ///< x | y | z: the mode travels along + this axis
	double position = 0;                ///< [m] injection plane along the propagation axis, from the structure centre
	std::string mpb_x_axis;             ///< FDTD axis of MPB's x, "": z for propagation x, x for y and z
	std::string mpb_y_axis;             ///< FDTD axis of MPB's y, "": y for x, z for y, y for z
	double offset_x = 0, offset_y = 0, offset_z = 0; ///< [m] MPB cell centre from the structure centre
	bool operator==(const ModeSpec&) const = default;
};
} // namespace fdtd::sources

#endif
