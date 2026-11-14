/** \file   materials/codes.h
 *  \brief  Material codes: how a code in the material maps selects a model and its parameter slot.
 *  @verbatim
      code < DispMateOffsetD (100)                       dielectric, slot = code
      100 .. 149   Drude + 2 critical points             slot mm = code - 100, mm < DrudeCP
      150 .. 199   Drude                                 mm < Drude
      200 .. 249   Lorentz                               mm < Lorentz
      250 .. 299   modified Lorentz                      mm < ModLorentz
      300 .. 349   hydrodynamic, nonlinear free electrons, mm < NonLinFreeElectron
      350 .. 399   Drude, hydrodynamic or two-temperature electron gas, + 1 Lorentz pole in the
                   critical-point arrays; mm < NonlinearDrudeLorentz, models/drude_lorentz.h
      400 .. 449   magnetic Lorentz                      mm < MagneticLorentz
      code + nl_code                                     a nonlinear cell of material `code`, SHG/THG
   @endverbatim
 */
#ifndef FDTD_MATERIALS_CODES_H
#define FDTD_MATERIALS_CODES_H

namespace fdtd::materials
{
inline constexpr int NRMAXMATERIALS = 600;	  ///< material codes 0 .. NRMAXMATERIALS - 1
inline constexpr int NrMaxCP = 4;			  ///< critical-point poles per material, at most
inline constexpr int nr_disp_mat_d = 500;	  ///< dispersive material slots
inline constexpr int DispMateOffsetD = NRMAXMATERIALS - nr_disp_mat_d; ///< first dispersive code (100)
// each constant is the END of its range, the start of the next: the Drude-Lorentz slots are
// NonLinFreeElectron .. NonlinearDrudeLorentz - 1, which is why parameters.cc offsets them by
// "DrudeLorentzShifter = NonLinFreeElectron"
inline constexpr int DrudeCP = 50;			  ///< end of the Drude + 2CP slots
inline constexpr int Drude = 100;			  ///< end of the Drude slots
inline constexpr int Lorentz = 150;			  ///< end of the Lorentz slots
inline constexpr int ModLorentz = 200;		  ///< end of the modified Lorentz slots
inline constexpr int NonLinFreeElectron = 250;	  ///< end of the hydrodynamic slots
inline constexpr int NonlinearDrudeLorentz = 300; ///< end of the nonlinear Drude-Lorentz slots
inline constexpr int MagneticLorentz = 350;		  ///< end of the magnetic Lorentz slots
inline constexpr int nl_code = 620;				  ///< added to the code of a nonlinear cell
} // namespace fdtd::materials

#endif // FDTD_MATERIALS_CODES_H
