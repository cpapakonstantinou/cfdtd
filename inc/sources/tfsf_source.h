/**
 * @file sources/tfsf_source.h
 * @brief The TF/SF plane-wave source of this rank: the injection into the 1D grid and the face list
 *
 * Built once per sweep from the run's flags. Each step:
 * @verbatim
      inject_h:     signal(t) into H of the 1D grid at the source cell, both lines with periodic x
                    and z, then the H corrections of this rank's TF/SF faces
      inject_e:     the E corrections of this rank's TF/SF faces
   @endverbatim
 * The 1D grid itself advances inside the E and H updates, sources/aux_grid.h.
 * Faces: the y faces always; the x faces unless x is periodic, H with pol z, E with pol x; the z
 * faces unless z is periodic, H with pol x, E with pol z.
 */
#ifndef FDTD_SOURCES_TFSF_SOURCE_H
#define FDTD_SOURCES_TFSF_SOURCE_H

#include "sources/signal.h"
#include "sources/source.h"
#include "sources/tfsf.h"
#include <vector>

namespace fdtd::sources
{
class TfsfSource : public Source
{
public:
	/// the plane wave of this rank, none when the run has no TF/SF plane wave, or uses the mode source
	static TfsfSource for_this_rank(const Signal &signal);

	bool active() const { return _active; }
	/// signal at time.h into the 1D grid, then the H face corrections
	void inject_h(const StepTime &time) const override;
	/// the E face corrections
	void inject_e(const StepTime &time) const override;

private:
	bool _active = false;
	std::vector<Face> _h_faces, _e_faces; ///< this rank's faces, in injection order
	const Signal *_signal = nullptr;      ///< the run's signal
};
} // namespace fdtd::sources

#endif
