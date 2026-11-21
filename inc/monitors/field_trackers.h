/**
 * @file monitors/field_trackers.h
 * @brief The field trackers: running DFT of E and H at the points of the multipole sphere
 *
 * The points are the cells at the integer positions
 * @verbatim
      int(R sin(theta_i) cos(phi_j)), int(R sin(theta_i) sin(phi_j)), int(R cos(theta_i))
      theta_i = pi i / N_theta,  phi_j = 2 pi j / N_phi, [multipoles] sphere_radius, sphere_theta, sphere_phi
   @endverbatim
 * around the shape centre; a point on the sphere grid that several (i, j) map to is tracked once per
 * (i, j). Each rank keeps the points it holds, in global (x, y, z) scan order.
 * @verbatim
      dft[set][E/H][component][point][frequency][cos, sin]
          set 0: the fields at analysis_omega[m];  set 1, separate linear/nonlinear fields: the harmonic fields at 2 analysis_omega[m]
   @endverbatim
 * Written at every checkpoint to frames/tracker_<n>.txt, and tracker_harmonic_<n>.txt, n the
 * sphere-grid index of the point.
 */
#ifndef FDTD_MONITORS_FIELD_TRACKERS_H
#define FDTD_MONITORS_FIELD_TRACKERS_H

#include "core/array.h"
#include <vector>

namespace fdtd::monitors
{
struct TrackPoint
{
	int x; ///< local index, with the inner origin
	int y; ///< local index, with the inner origin
	int z; ///< local index, with the inner origin
	int index;   ///< index on the sphere grid (N_theta x N_phi)
};

struct FieldTrackers
{
	std::vector<TrackPoint> point;
	array3<double> dft[2][2][3];

	/// finds this rank's points around the centre and allocates their DFTs, collective: prints the total
	void setup(int center_x, int center_y, int center_z);
	/// accumulates this step: E at t_e, H at t_h
	void sample(double t_e, double t_h);
	/// writes one file per point, and per harmonic set
	void write() const;
	/// forgets the points and frees their DFTs
	void release();
};

FieldTrackers &field_trackers();
} // namespace fdtd::monitors

#endif
