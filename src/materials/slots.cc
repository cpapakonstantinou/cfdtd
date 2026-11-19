/**
 * @file materials/slots.cc
 * @brief The run's material slots, materials/slots.h: each definition into the parameter tables
 */
#include "materials/slots.h"

#include <cmath>
#include <stdexcept>

#include "core/constants.h"
#include "materials/dispersive_state.h"
#include "materials/hydro/state.h"
#include "materials/nonlinear/state.h"
#include "materials/thermal/thermal.h"

namespace fdtd::materials
{
Slots the_slots;

int first_slot(Kind k)
{
	for (int s = 0; s < the_slots.count; ++s)
		if (the_slots.kind[s] == k)
			return s;
	return -1;
}

namespace
{
/// the nonlocal pressure coefficient of the hydrodynamic model in 3D, beta = v_F sqrt(C / D)
double nonlocal_coefficient()
{
	const double d = 3.0;
	const double c = 3.0 * d / (d + 2.0);
	return std::sqrt(c / d);
}

/// the poles of a critical-point or Lorentz-pole model, phi_cp holds the phase or the numerator term
void load_poles(ModelParameters &mp, int s, const Definition &d)
{
	mp.nr_cp[s] = double(d.poles.size());
	for (size_t p = 0; p < d.poles.size(); ++p)
	{
		const Pole &q = d.poles[p];
		mp.big_omega_cp[s][p] = q.frequency;
		mp.big_gamma_cp[s][p] = q.damping;
		mp.a_coeff_cp[s][p] = q.amplitude;
		mp.phi_cp[s][p] = d.kind == Kind::modified_lorentz ? q.numerator : q.phase;
	}
}

/// the free-electron parameters every electron-gas model reads
void load_electron_gas(ModelParameters &mp, int s, const Definition &d, const LoadSpec &spec)
{
	mp.beta_nl[s] = d.fermi_velocity * spec.fermi_velocity_scale * nonlocal_coefficient();
	mp.diff_coeff[s] = d.diffusion;
	mp.equilibrium_electron_density[s] = d.electron_density;
}

/// chi2 or chi3 of slot s, and whether the slot is nonlinear
void load_tensors(int s, const Definition &d, const LoadSpec &spec)
{
	namespace nl = fdtd::materials::nonlinear;
	for (int i = 0; i < 18; ++i)
		nl::chi_2[s][i] = spec.harmonic_order == 2 ? double(2) * d.chi2_values[i] : 0; // chi2 = 2 d
	for (int i = 0; i < 21; ++i)
		nl::chi_3[s][i] = spec.harmonic_order == 3 ? d.chi3_values[i] : 0;
	auto diagonal = [](const double *t, int n) {
		bool any = false;
		for (int i = 0; i < n; i += 7)
			any = any || t[i] != 0;
		return any;
	};
	auto off_diagonal = [](const double *t, int n) {
		bool any = false;
		for (int i = 0; i < n; ++i)
			if (i % 7 != 0)
				any = any || t[i] != 0;
		return any;
	};
	const bool d2 = diagonal(nl::chi_2[s], 18);
	const bool o2 = off_diagonal(nl::chi_2[s], 18);
	const bool d3 = diagonal(nl::chi_3[s], 21);
	const bool o3 = off_diagonal(nl::chi_3[s], 21);
	nl::nl[s] = d2 || o2 || d3 || o3;
	nl::iso_chi_2[s] = d2 && !o2;
	nl::iso_chi_3[s] = d3 && !o3;
}
} // namespace

void load_slots(const std::vector<Definition> &definitions, const LoadSpec &spec)
{
	using namespace fdtd::constants;
	if (definitions.empty() || definitions[0].name != "vacuum")
		throw std::logic_error("load_slots: slot 0 is vacuum");
	if (int(definitions.size()) > max_slots)
		throw std::invalid_argument("the run uses " + std::to_string(definitions.size()) + " materials, at most " + std::to_string(max_slots));
	ModelParameters &mp = model_parameters();
	mp = ModelParameters{};
	the_slots = Slots{};
	the_slots.count = int(definitions.size());
	bool add_term_set = false;
	for (int s = 0; s < the_slots.count; ++s)
	{
		const Definition &d = definitions[s];
		the_slots.kind[s] = d.kind;
		the_slots.names.push_back(d.name);
		mp.sig_m[s] = d.conductivity;
		load_tensors(s, d, spec);
		switch (d.kind)
		{
		case Kind::dielectric:
			mp.epre[s] = d.n != 0 ? d.n * d.n : d.eps;
			break;
		case Kind::drude_cp:
		{
			mp.epre_inf[s] = d.eps_inf;
			mp.gam_plasma[s] = d.plasma_damping;
			if (d.sheet_thickness > 0)
			{
				const double f = d.sheet_thickness / spec.dy;
				mp.ome_plasma[s] = d.plasma_frequency * std::sqrt(d.sheet_scale * f);
				load_poles(mp, s, d);
				for (size_t p = 0; p < d.poles.size(); ++p)
					mp.a_coeff_cp[s][p] = d.poles[p].amplitude * d.sheet_scale * f;
			}
			else
			{
				mp.ome_plasma[s] = d.plasma_frequency;
				load_poles(mp, s, d);
			}
			load_electron_gas(mp, s, d, spec);
			break;
		}
		case Kind::drude:
			mp.epre_inf[s] = d.eps_inf;
			mp.gam_plasma[s] = d.plasma_damping;
			if (d.carrier_density > 0)
				mp.ome_plasma[s] = std::sqrt((d.carrier_density * E0 * E0) / (EP0 * d.effective_mass * ME));
			else
				mp.ome_plasma[s] = d.plasma_frequency;
			load_electron_gas(mp, s, d, spec);
			break;
		case Kind::lorentz:
		case Kind::magnetic_lorentz:
			mp.epre_sta[s] = d.eps_static;
			mp.epre_inf[s] = d.eps_inf;
			mp.ome_plasma[s] = d.resonance_frequency;
			mp.gam_plasma[s] = d.damping;
			mp.a_3[s] = d.magnetic_coupling;
			break;
		case Kind::modified_lorentz:
			mp.epre_inf[s] = d.eps_inf;
			load_poles(mp, s, d);
			break;
		case Kind::hydrodynamic:
			mp.epre_inf[s] = d.eps_inf;
			mp.ome_plasma[s] = d.plasma_frequency;
			mp.gam_plasma[s] = d.plasma_damping;
			load_poles(mp, s, d);
			load_electron_gas(mp, s, d, spec);
			mp.effective_electron_mass[s] = d.effective_mass > 0 ? d.effective_mass * ME
															  : E0 * E0 * d.electron_density / (mp.ome_plasma[s] * mp.ome_plasma[s] * EP0);
			break;
		case Kind::drude_lorentz:
			mp.epre_inf[s] = d.eps_inf;
			mp.ome_plasma[s] = d.plasma_frequency;
			mp.gam_plasma[s] = d.plasma_damping;
			load_poles(mp, s, d);
			load_electron_gas(mp, s, d, spec);
			if (!add_term_set)
			{
				// the two-temperature plasma frequency of the first Drude-Lorentz material starts at its own
				mp.plasma_freq_add_term = 0;
				mp.plasma_freq_add_term = mp.ome_plasma[s] - fdtd::materials::thermal::plasma_frequency(DEFAULT_TEMPERATURE);
				add_term_set = true;
			}
			break;
		}
	}
}
} // namespace fdtd::materials
