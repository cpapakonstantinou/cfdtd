/**
 * @file materials/permittivity.cc
 * @brief The relative permittivity and refractive index of a slot at a frequency, materials/materials.h
 *
 * The two differ in the Lorentz model: the refractive index takes 2 i omega gamma in the denominator,
 * the permittivity i omega gamma, as the solver always had. docs/dev/Roadmap.md §7 records it.
 */
#include "materials/materials.h"

#include <cmath>

#include "core/constants.h"
#include "materials/dispersive_state.h"
#include "materials/slots.h"

namespace fdtd::materials
{
namespace
{
/// eps_r of slot s at omega; `lorentz_factor` multiplies the damping term of the Lorentz model
std::complex<double> permittivity_of(int s, double omega, double lorentz_factor)
{
	using fdtd::constants::I_UNIT;
	const auto &mp = model_parameters();
	std::complex<double> eps = std::complex<double>(0.);
	switch (kind(s))
	{
	case Kind::dielectric:
		eps = std::complex<double>(mp.epre[s], 0.);
		break;
	case Kind::drude_cp:
		eps = mp.epre_inf[s] - mp.ome_plasma[s] * mp.ome_plasma[s] / omega / (omega + I_UNIT * mp.gam_plasma[s]);
		for (int p = 0; p < mp.nr_cp[s]; ++p)
			eps += mp.a_coeff_cp[s][p] * mp.big_omega_cp[s][p] *
				   (std::complex<double>(exp(std::complex<double>(+I_UNIT * mp.phi_cp[s][p]))) / (mp.big_omega_cp[s][p] - omega - I_UNIT * mp.big_gamma_cp[s][p]) +
					std::complex<double>(exp(std::complex<double>(-I_UNIT * mp.phi_cp[s][p]))) / (mp.big_omega_cp[s][p] + omega + I_UNIT * mp.big_gamma_cp[s][p]));
		break;
	case Kind::drude:
		eps = mp.epre_inf[s] - mp.ome_plasma[s] * mp.ome_plasma[s] / omega / (omega + I_UNIT * mp.gam_plasma[s]);
		break;
	case Kind::lorentz:
		eps = mp.ome_plasma[s] * mp.ome_plasma[s] / (mp.ome_plasma[s] * mp.ome_plasma[s] - lorentz_factor * I_UNIT * omega * mp.gam_plasma[s] - omega * omega);
		eps *= (mp.epre_sta[s] - mp.epre_inf[s]);
		eps += mp.epre_inf[s];
		break;
	case Kind::modified_lorentz:
		eps = mp.epre_inf[s];
		for (int p = 0; p < mp.nr_cp[s]; ++p)
			eps += mp.a_coeff_cp[s][p] * (mp.big_omega_cp[s][p] * mp.big_omega_cp[s][p] - I_UNIT * omega * mp.phi_cp[s][p]) /
				   (mp.big_omega_cp[s][p] * mp.big_omega_cp[s][p] - double(2.) * I_UNIT * omega * mp.big_gamma_cp[s][p] - omega * omega);
		break;
	default: // the models without a closed form here read 0
		break;
	}
	return eps;
}
} // namespace

std::complex<double> refractive_index(int s, double omega)
{
	const std::complex<double> eps = permittivity_of(s, omega, double(2.));
	const double r = sqrt(eps.real() * eps.real() + eps.imag() * eps.imag());
	return std::complex<double>(sqrt((r + eps.real()) / 2), sqrt((r - eps.real()) / 2));
}

std::complex<double> complex_permittivity(int s, double omega) { return permittivity_of(s, omega, double(1.)); }
} // namespace fdtd::materials
