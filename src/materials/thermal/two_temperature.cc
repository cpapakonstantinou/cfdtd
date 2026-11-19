#include "grid/subdomain.h"
#include "materials/state.h"
#include "fields/fields.h"
#include "materials/hydro/state.h"
#include "io/text.h"
#include "materials/thermal/state.h"
#include "materials/thermal/thermal.h"
#include "grid/halo_exchange.h"
#include "materials/spec.h"
using fdtd::fields::Ex;
using fdtd::fields::Ey;
using fdtd::fields::Ez;
using fdtd::grid::subdomain;

namespace
{
} // namespace
using namespace fdtd::materials::thermal; // transitional: module state, thermal/state.h

void fdtd::materials::thermal::update_ttm_rk(){
  auto &ds = fdtd::materials::dispersive_state();
  int i;
  int j;
  int k;
  int ii;
  int jj;
  int kk;
  double Ce;
  double EP_Coup;
  double Cl;
  double Te;
  double Tl;
  double P;
  int nrk;
  double y_vec[4][3];
  double Func[4][3];

  for (i = 0; i <= subdomain.cells_x; ++i) {
    ii=i+1;

    for (j = 0; j <= subdomain.cells_y; ++j) {
      jj=j+1;

      for (k = 0; k <= subdomain.cells_z; ++k) {
        kk=k+1;

  Te = electron_temp[i][j][k];
  Tl = lattice_temp[i][j][k];

  y_vec[0][0] = Te;
  y_vec[0][1] = Tl;

  P = ((fdtd::materials::hydro::ex_2[ii][jj][kk])*(ds.c[0].j[fdtd::materials::thermal::heated_slot()][ii][jj][kk]) + (fdtd::materials::hydro::ey_2[ii][jj][kk])*(ds.c[1].j[fdtd::materials::thermal::heated_slot()][ii][jj][kk]) + (fdtd::materials::hydro::ez_2[ii][jj][kk])*(ds.c[2].j[fdtd::materials::thermal::heated_slot()][ii][jj][kk]));

    for(nrk = 0; nrk < num_rk; nrk++){

      Ce = fdtd::materials::thermal::electron_heat_capacity(y_vec[0][0]);
      Cl = fdtd::materials::thermal::lattice_heat_capacity(y_vec[0][1]);
      EP_Coup = fdtd::materials::thermal::electron_phonon_coupling(y_vec[0][0]);

      Func[0][0] = fdtd::materials::model_spec().dt/num_rk * (-EP_Coup*(y_vec[0][0] - y_vec[0][1]) + P)/Ce;
      Func[0][1] = fdtd::materials::model_spec().dt/num_rk * (EP_Coup*(y_vec[0][0] - y_vec[0][1]))/Cl;

      y_vec[1][0] = y_vec[0][0] + 0.5*Func[0][0];
      y_vec[1][1] = y_vec[0][1] + 0.5*Func[0][1];

      Ce = fdtd::materials::thermal::electron_heat_capacity(y_vec[1][0]);
      Cl = fdtd::materials::thermal::lattice_heat_capacity(y_vec[1][1]);
      EP_Coup = fdtd::materials::thermal::electron_phonon_coupling(y_vec[1][0]);

      Func[1][0] = fdtd::materials::model_spec().dt/num_rk * (-EP_Coup*(y_vec[1][0] - y_vec[1][1]) + P)/Ce;
      Func[1][1] = fdtd::materials::model_spec().dt/num_rk * (EP_Coup*(y_vec[1][0] - y_vec[1][1]))/Cl;

      y_vec[2][0] = y_vec[0][0] + 0.5*Func[1][0];
      y_vec[2][1] = y_vec[0][1] + 0.5*Func[1][1];

      Ce = fdtd::materials::thermal::electron_heat_capacity(y_vec[2][0]);
      Cl = fdtd::materials::thermal::lattice_heat_capacity(y_vec[2][1]);
      EP_Coup = fdtd::materials::thermal::electron_phonon_coupling(y_vec[2][0]);

      Func[2][0] = fdtd::materials::model_spec().dt/num_rk * (-EP_Coup*(y_vec[2][0] - y_vec[2][1]) + P)/Ce;
      Func[2][1] = fdtd::materials::model_spec().dt/num_rk * (EP_Coup*(y_vec[2][0] - y_vec[2][1]))/Cl;

      y_vec[3][0] = y_vec[0][0] + 0.5*Func[2][0];
      y_vec[3][1] = y_vec[0][1] + 0.5*Func[2][1];

      Ce = fdtd::materials::thermal::electron_heat_capacity(y_vec[3][0]);
      Cl = fdtd::materials::thermal::lattice_heat_capacity(y_vec[3][1]);
      EP_Coup = fdtd::materials::thermal::electron_phonon_coupling(y_vec[3][0]);

      Func[3][0] = fdtd::materials::model_spec().dt/num_rk * (-EP_Coup*(y_vec[3][0] - y_vec[3][1]) + P)/Ce;
      Func[3][1] = fdtd::materials::model_spec().dt/num_rk * (EP_Coup*(y_vec[3][0] - y_vec[3][1]))/Cl;

      y_vec[0][0] += (Func[0][0] + 2.0*Func[1][0] + 2.0*Func[2][0] + Func[3][0])/6.0;
      y_vec[0][1] += (Func[0][1] + 2.0*Func[1][1] + 2.0*Func[2][1] + Func[3][1])/6.0;

    }

    electron_temp[i][j][k] = y_vec[0][0];
    lattice_temp[i][j][k] = y_vec[0][1];

      }
    }
  }
}

bool fdtd::materials::thermal::update_ttm(int step)
{
    auto &ds = fdtd::materials::dispersive_state();
    bool OUTPUT_PUMP_DATA2 = false;
    int i;
    int j;
    int k;
    int ii;
    int jj;
    int kk;
    for (i = 0; i <= subdomain.cells_x; i++) {
        ii = i + 1;
        for (j = 0; j <= subdomain.cells_y; j++) {
            jj = j + 1;
            for (k = 0; k <= subdomain.cells_z; k++) {
                kk = k + 1;
                fdtd::materials::hydro::ex_2[ii][jj][kk] = Ex[i][j][k];
                fdtd::materials::hydro::ey_2[ii][jj][kk] = Ey[i][j][k];
                fdtd::materials::hydro::ez_2[ii][jj][kk] = Ez[i][j][k];
            }
        }
    }

    fdtd::grid::exchange_margins(fdtd::materials::hydro::ex_2, fdtd::materials::hydro::ey_2, fdtd::materials::hydro::ez_2);
    fdtd::grid::exchange_margins(ds.c[0].j[heated_slot()], ds.c[1].j[heated_slot()], ds.c[2].j[heated_slot()]); // only for ITO

    if (!fdtd::materials::model_spec().probe)
        fdtd::materials::thermal::update_ttm_rk();

    if (fdtd::materials::model_spec().pump_probe && fdtd::materials::model_spec().pump) {
        i = int(subdomain.cells_x / 2);
        k = int(subdomain.cells_z / 2);
        for (j = 0; j <= subdomain.cells_y; j++) {
            temp_save[0][step][j] = electron_temp[i][j][k];
            if (electron_temp[i][j][k] > 301.0) {
                OUTPUT_PUMP_DATA2 = true;

            }
        }
    }

    if (fdtd::materials::model_spec().pump_probe && fdtd::materials::model_spec().probe) {
        for (i = 0; i <= subdomain.cells_x; i++) {
            for (j = 0; j <= subdomain.cells_y; j++) {
                for (k = 0; k <= subdomain.cells_z; k++) {
                    if (pump_probe_delay < 0) {
                        if (step < -int(pump_probe_delay / fdtd::materials::model_spec().dt))
                            electron_temp[i][j][k] = 300.0;
                        else
                            electron_temp[i][j][k] = temp_save[0][step + int(pump_probe_delay / fdtd::materials::model_spec().dt)][j];
                    }
                    else {
                        if (step + int(pump_probe_delay / fdtd::materials::model_spec().dt) < nr_max_time_steps_pump)
                            electron_temp[i][j][k] = temp_save[0][step + int(pump_probe_delay / fdtd::materials::model_spec().dt)][j];
                        else
                            electron_temp[i][j][k] = electron_temp[i][j][k];
                    }
                }
            }
        }
    }
    return OUTPUT_PUMP_DATA2;

}

double fdtd::materials::thermal::electron_phonon_coupling(double Te){
  return (-2.25621369e-05*Te*Te*Te*Te*Te + 3.50483834e+00*Te*Te*Te*Te -2.12068062e+05*Te*Te*Te + 6.38479102e+09*Te*Te + 2.17060773e+13*Te + 9.00926467e+17)/7.0;

}

double fdtd::materials::thermal::electron_heat_capacity(double Te){
  double Tesqrt = sqrt(Te);

  return  (-9.79573687e-07*Tesqrt*Tesqrt*Tesqrt*Tesqrt*Tesqrt + 6.30170887e-04*Tesqrt*Tesqrt*Tesqrt*Tesqrt - 1.39337340e-01*Tesqrt*Tesqrt*Tesqrt + 9.87737302e+00*Tesqrt*Tesqrt + 4.75208334e+02*Tesqrt - 7.23459681e+03);
}

double fdtd::materials::thermal::lattice_heat_capacity([[maybe_unused]] double Tl){
  return 2.54e6;
}

double fdtd::materials::thermal::plasma_frequency(double Te){
  const auto &mp = fdtd::materials::model_parameters();
  return 2.16439846e-08*Te*Te*Te*Te*Te - 2.76226816e-03*Te*Te*Te*Te + 1.12417870e+02*Te*Te*Te - 6.94572643e+05*Te*Te -6.94024494e+10*Te + 2.71693490e+15 + mp.plasma_freq_add_term;

}

double fdtd::materials::thermal::damping_constant(double Te, [[maybe_unused]] double Tl){
	const auto &mp = fdtd::materials::model_parameters();

  return mp.gam_plasma[heated_slot()]*0.98*(1.0 + (Te/15000.0));

}
