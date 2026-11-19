# Epsilon Near Zero, Nonlinear optical modelling of Transparent Conducting Oxides

Transparent conducting oxides (TCOs) are optically transparent materials with large electrical conductivities. They are wide-bandgap, degenerate semiconductors that possess a metallic optical response in the near infrared (NIR) due to an abundance of free electrons in the conduction band. Thus their optical properties, in the NIR are best described using the Drude model

$$
	\frac{\partial\mathbf{J}}{\partial t}=-\gamma\mathbf{J}+\varepsilon_0\omega_p^2\mathbf{E}
$$

 where $\omega_p=\sqrt{Ne^2/\varepsilon_0m^\ast}$ is the plasma frequency, $\gamma$ is the damping coefficient, $m^{*}$ is the effective mass of the conduction electrons, $e$ is the charge of an electron, $N$ is the free electron density, and $\varepsilon_0$ is the permittivity of free space. A static permittivity $\varepsilon_{\infty}$ or a full Lorentz model can be used for the bound electrons. For most metals, the conduction band is mostly parabolic and thus the $m^{*}$  (and $\omega_p$) are constant. This is generally not the case in TCOs exposed to intense, ultrafast laser pulses due to their high doping concentration (and thus high carrier density). The conduction electrons are pushed to the non-parabolic region of the conduction band. As a result $m^{*}$ and subsequently $\omega_p$ becomes dependent on the electron energies, and thus the electron temperature $T_e(\textbf{r},t)$. This is why we call TCOs an active/dynamic nanophotonic material; it's plasma frequency is not constant, but will vary in time and space when irradiated with an intense, ultrafast pulse. The functional form of the $\omega_p$

 $$    
     \omega_p^2\left(T_e\right)=\frac{e^2}{3m_0^\ast\varepsilon_0\pi^2}\left(\frac{2m_0^\ast}{\hbar^2}\right)^{3/2}\int_{0}^{\infty}{dE\ {(E+CE^2)}^{3/2}\left(1+2CE\right)^{-1}\left(-\frac{\partial f_{FD}\left(E,T_e\right)}{\partial E}\right)}
$$

where $E$ is the electron energy, $\hbar$ is the reduced Planck’s constant, $m_0^*$ is the effective mass at the conduction band minimum, and $C$ is a non-parabolicity factor, $f_{FD}=(\exp{((E-\mu)/k_B T_e)} + 1)^{-1}$ is the Fermi-Dirac distribution, $\mu$ is the chemical potential, and $k_B$ is the Boltzmann constant.

Our code includes a model for the most commonly used TCO, indium tin oxide (ITO). During the FDTD simulation, where Maxwell's equations are solved along with the Drude model,

$$    
	\mu_0\frac{\partial\mathbf{H}}{\partial t}=-\mathrm{\nabla}\times\mathbf{E} \label{HM},\\
	{\varepsilon_\infty\varepsilon}_0\frac{\partial\mathbf{E}}{\partial t}=\nabla\times\mathbf{H}-\mathbf{J} \label{EM},\\
	\frac{\partial\mathbf{J}}{\partial t}=-\gamma\mathbf{J}+\varepsilon_0\omega_p^2(T_e)\mathbf{E}, \label{Current}
$$

where $T_e$ is calculated by numerical calculation of the two-temperature model (TTM)

$$    
	C_e\left(T_e\right)\frac{\partial T_e(t,\mathbf{r})}{\partial t}&=-g_{ep}\left(T_e,T_l\right)\left(T_e\left(t,\mathbf{r}\right)-T_l\left(t,\mathbf{r}\right)\right)\ +P\left(\mathbf{r},t\right),\\
	C_l\left(T_l\right)\frac{\partial T_l\left(t,\mathbf{r}\right)}{\partial t}&=g_{ep}\left(T_e,T_l\right)\left(T_e\left(t,\mathbf{r}\right)-T_l\left(t,\mathbf{r}\right)\right).
$$

where $T_l$ is the lattice temperature, $C_e$ is the electron heat capacity, $C_l$ is the lattice heat capacity, $g_{ep}$ is the electron-phonon coupling coefficient, and $P$ is the electromagnetic power absorbed by the conduction electrons.

The TTM and electrodynamics equations are coupled via two terms:  (i) the absorbed power $P=\textbf{J}\cdot\textbf{E}$ given by Poynting’s theorem, and used as the source in TTM and (ii) the plasma frequency $\omega_p(T_e)$ is used as the source of the nonlinearity in the electrodynamics equations. Maxwell's Equations + Drude model are calculated via FDTD, and the TTM is updated at each FDTD step via a 4th order Runge-Kutta scheme.

For derivations, numerical values, and more details, please see our pre-print [1] and the thesis of Luis Enrique Cortes Herrera [2].

## Pump-Probe simulations

Because this is a time-varying optical material, it is common to run pump-probe experiments to
determine the lifetime of the nonlinear processes. Both methods below need separate simulations for
the pump and the probe pulse.

1. **1D structures (thin films): a true pump-probe simulation.** The pump run
   (`[temperature] model = pump`) saves the temperature field at each time step and each position
   along the propagation direction y. The probe run (`model = probe`) loads it: the probe pulse
   has no effect on the temperature, it only feels the temperature left by the pump.
   `[pump_probe] delay` [s] is the delay between the two pulses. One pump run serves any number of
   probe runs.

2. **Higher-dimensional structures: a frozen temperature.** The memory needed to keep $T_e$ for all
   times and positions is too large. Instead, run the pump with `model = two_temperature` and
   `[two_temperature] save_frames = true`; the 3D temperature field is written to
   `state/electron_temperature_<step>.h5` every N steps, N being the largest frame period of
   `[output] frames_xy`, `frames_xz`, `frames_yz`, `frames_3d` (N = 10000 is typical for
   dx = 2 nm). For the probe run, copy the file of the time step to probe to
   `<input_dir>/electron_temperature.h5` and set `model = load`: that temperature field stays for
   the whole run. This probes the optical properties at exactly that time step, so it is not a
   true pump-probe calculation like method 1.

## How to use

Select a `drude_lorentz` material, [Materials](Materials.md); the two-temperature model heats the
first one of the run. The optical properties of thin-film ITO depend on the sample and its
thickness, and the catalog has one material per film: `ito_film_thin` the parameters of [3],
reported for 30 nm, `ito_film_mid` those of [4], reported for 220 nm, and `ito_film_thick` those of
[5], reported for 310 nm. A film of another sample is a definition of its own.

```ini
[geometry]
middle_medium = ito_film_thick
middle_thickness = 310e-9

[temperature]
model = two_temperature

[two_temperature]
rk_stages = 1            ; Runge-Kutta stages of the temperature update
save_frames = false      ; write temperature frames
```

| `[temperature] model` | Run |
|---|---|
| `none` (default) | no two-temperature model: the plain Drude-Lorentz response |
| `two_temperature` | Drude-Lorentz + two-temperature model with $\omega_p(T_e)$ |
| `load` | probe with a frozen temperature field from `<input_dir>/electron_temperature.h5` (method 2) |
| `pump` | pump run of a 1D pump-probe pair (method 1) |
| `probe` | probe run of a 1D pump-probe pair; `[pump_probe] delay` |

The two-temperature model combines with `[electrons] model = local` or `hydrodynamic`. A run with
it cannot be continued from a state file yet ([Stopping and Restarting](../tutorial/Restart.md)).

## References

[1] doi.org/10.48550/arXiv.2110.14806

[2] Ask Prof. Israel De Leon Arizpe, or Prof. Lora Ramunno for this

[3] doi.org/10.1038/s41566-017-0089-9

[4] doi.org/10.1103/PhysRevApplied.11.064062

[5]	doi.org/10.1126/science.aae0330
