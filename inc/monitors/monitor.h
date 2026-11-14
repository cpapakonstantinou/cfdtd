/**
 * @file monitors/monitor.h
 * @brief The monitors of a run as one list of objects with a sampling hook
 *
 * After each full step the time loop calls sample on every monitor in list order, continuous-wave
 * runs only after their settling periods. The list is built once per sweep (build_monitors) from
 * the run's configuration, so no flag is tested inside the loop.
 * @verbatim
      monitor           samples                                  finish, output,               when
      volume DFT        E, and H, over the inner domain          exact multipoles; at the end  field saves, absorption,
                                                                 3D E/H, field saves
                                                                                               DFG, exact multipoles
      P_NL DFT          nonlinear polarisation, harmonic band    -                             exact multipoles, SHG/THG
      ART planes        E, H on the reflection/transmission      R, T, A;           also at    always, planes where present
                        planes                                   steps N - kN/10, k = 1..9
      near-to-far       -                                        scattering-box faces (HDF5)   SHG/THG with two frequencies
      flux boxes        E, H on the faces of the Q boxes         Q coefficients, ART-TFSF      q_coeff
                                                                 also at checkpoints
      field points      -                                        near/far-field point files    always
      display planes    E (H) on the XY / YZ / XZ planes         plane maps (HDF5)             always, planes where present
      field trackers    E, H at the tracking points              -                             multipole sphere, MULTIPOLES
      field trackers, checkpoints: one file per point
      custom monitors   DFT boxes and time probes from the INI   at the end: DFT boxes         [monitor] sections
                                                                 probes written in sample
      time probes       E, H at single cells, every step         checkpoints: trace files      [probe] sections
   @endverbatim
 * E is sampled at the E time of the step, H at the H time, half a step earlier. After the loop the
 * driver calls normalize on every monitor, runs with a DFT, not dipole runs, before anything reads the
 * spectra, multipoles, ART, near-to-far, Q coefficients, field monitors, plots.
 */
#ifndef FDTD_MONITORS_MONITOR_H
#define FDTD_MONITORS_MONITOR_H

#include <memory>
#include <vector>

namespace fdtd::monitors
{
/// the times a monitor samples at in one step
struct SampleTime
{
	double e;  ///< time of E
	double h;  ///< time of H
	int step;  ///< time step, probes
};

class Monitor
{
public:
	virtual ~Monitor() = default;
	virtual void sample(const SampleTime &time) = 0;
	/// sampled also during the settling periods of continuous-wave runs, time-domain monitors
	virtual bool every_step() const { return false; }
	/// after the time loop: divide the running DFTs by the source spectrum, nothing for time-domain monitors
	virtual void normalize() {}
	/// after normalize: the monitor's analysis and output files. Often collective, so the list holds
	/// the same monitors, in the same order, on every rank, rank-local work is skipped inside sample
	virtual void finish() {}
	/// once at the end of the run, after the last sweep, or instead of the time loop on the mode-source
	/// dummy run: the outputs kept per run, not per sweep; collective like finish
	virtual void write_final() {}
	/// after every step, before the checkpoint test: output on the monitor's own schedule; collective
	virtual void after_step(int step) { (void)step; }
	/// at every checkpoint, [run] checkpoint_count per run: intermediate output; collective like finish
	virtual void checkpoint(int step) { (void)step; }
};

/// the monitors of the run, in sampling order, the same list on every rank
std::vector<std::unique_ptr<Monitor>> build_monitors();
} // namespace fdtd::monitors

#endif
