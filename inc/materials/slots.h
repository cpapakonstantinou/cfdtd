/**
 * @file materials/slots.h
 * @brief The run's material slots: the materials in use, numbered from 0, and what each one is
 *
 * The material maps hold a slot number per point, the parameter and coefficient tables and the
 * per-material state are indexed by it. Slot 0 is vacuum. A run loads its definitions once,
 * load_slots; the model of a slot is kind(s).
 */
#ifndef FDTD_MATERIALS_SLOTS_H
#define FDTD_MATERIALS_SLOTS_H

#include <string>
#include <vector>

#include "materials/definition.h"

namespace fdtd::materials
{
inline constexpr int max_slots = 64; ///< materials in one run, at most

/// what the loader needs from the run besides the definitions
struct LoadSpec
{
	int harmonic_order = 0;          ///< 2 or 3: the tensors of that order are loaded, 0: linear run
	double dy = 0;                   ///< [m] cell size along y, sheet materials
	double fermi_velocity_scale = 1; ///< [electrons] fermi_velocity
};

/// the slot table: the model and the name of every slot
struct Slots
{
	int count = 0;
	Kind kind[max_slots] = {};
	std::vector<std::string> names;
};

extern Slots the_slots;

/// the model of slot s
inline Kind kind(int s) { return the_slots.kind[s]; }
/// slot s is a dispersive material
inline bool dispersive(int s) { return the_slots.kind[s] != Kind::dielectric; }
/// slot s holds model k
inline bool is(int s, Kind k) { return the_slots.kind[s] == k; }
/// the first slot of model k, -1 when none
int first_slot(Kind k);
/// slots in use
inline int slot_count() { return the_slots.count; }
inline const std::string &slot_name(int s) { return the_slots.names[s]; }

/// fills the slot table and the parameter tables, slot s from definitions[s]; definitions[0] is vacuum
void load_slots(const std::vector<Definition> &definitions, const LoadSpec &spec);
} // namespace fdtd::materials

#endif
