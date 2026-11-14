# Nonlocal Plasmonic Modelling

## Uses
This model is used for linear, nonlocal modelling of small plasmonic particles

## Model Details
Most optical models are considered "local" models as they relate the polarization field at a given point $\textbf{P(r)}$ to the local electric field at the same point $\textbf{E(r)}$. These models can be thought of as dipole approximations. Nonlocal models have higher-order multipole terms; $\textbf{P(r)}$ depends on the field divergence.

The most-used nonlocal plasmonics model is the nonlocal Drude model

$$
    \frac{\partial^2 \textbf{P}_{NL}}{\partial t^2} + \gamma \frac{\partial \textbf{P}_{NL}}{\partial t}=\varepsilon_0\omega_p^2\textbf{E}+\beta^2\nabla (\nabla\cdot\textbf{P}_{NL}),
$$

where $P_{NL}$ is the (nonlocal) free electron polarization field, $\omega_p$ is the plasma frequency, $\gamma$ is a collisional damping rate, and $\beta^2=1/3 v_F^2$ and $v_F$ is the Fermi velocity. What makes this model nonlocal is the last term which we call the "pressure term" as it describes the electron degeneracy pressure.

This model need only be used in the case of small particles and feature sizes. The pressure term is typically only non-zero near the boundaries of small particles (or sharp features) where there will be a large build up of charge in the presence of an electric field. This charge build-up only occurs within the Fermi wavelength $\lambda_F\approx 0.5$ nm of the boundary, and this the model requires a step size of $\Delta x < \lambda_F$.

As a further correction to the models, we can also account for charge density diffusion which changes how we define the current density (see Eq. 22 of [1] and results in the Generalized Nonlocal Optical Response (GNOR) model

$$
    \frac{\partial^2 \textbf{P}_G}{\partial t^2} + \gamma \frac{\partial \textbf{P}_G}{\partial t}=\varepsilon_0\omega_p^2\textbf{E}+\eta\nabla (\nabla\cdot\textbf{P}_G) + \mathcal{D}\frac{\partial}{\partial t}\nabla (\nabla\cdot\textbf{P}_G),
$$

where $P_{G}$ is the (GNOR) free electron polarization field, $\eta=\beta^2+\mathcal{D}\gamma$, and where $\mathcal{D}$ is the diffusion coefficient.

These models are programmed for Drude + 2 Critical Point materials, and as such the local 2CP model is used to accound for the bound electrons. For derivation, implementation details, and numerical values (which are hard-coded already), see our publication [2] and the review article [1].

## How to use

Select a Drude + critical points material (codes 100 to 149: 100 is gold, 101 silver; see
[Materials](Materials.md)) and the electron model:

```ini
[geometry]
material = 100

[electrons]
model = nonlocal         ; local (default) | nonlocal | gnor | hydrodynamic
fermi_velocity = 1.0     ; scales the Fermi velocity of the material
diffusion = 3.61e-4      ; diffusion coefficient [m^2/s], used by gnor
```

| `model` | Model |
|---|---|
| `local` | the local Drude + critical points response (default) |
| `nonlocal` | the nonlocal Drude model |
| `gnor` | the generalised nonlocal optical response: nonlocal with electron diffusion |
| `hydrodynamic` | the nonlinear hydrodynamic model ([Hydrodynamics](Hydrodynamics.md)) |

## References
[1] doi:10.1088/0953-8984/27/18/183204

[2] doi:10.1109/TAP.2020.3044579
