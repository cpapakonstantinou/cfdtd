# Hydrodynamic Plasma Model

## Uses


## The Model
The hydrodynamic plasma model is a nonlinear partial differential equation that describes how free electrons respond to electromagnetic fields. In it's simplest form, it is written as an equation of the macroscopic velocity field of the free electrons

$$
m\Bigg(\frac{\partial}{\partial t} + \textbf{v}\cdot\nabla + \gamma\Bigg)\textbf{v} = -e(\textbf{E}+\textbf{v}\times\textbf{B}) -\nabla\frac{\delta G[n]}{\delta n}
$$

where $\textbf{v}$ is the velocity field, $\textbf{E}$ and $\textbf{B}$ are the electric and magnetic fields, $m$ is the mass of an electron, $\gamma$ is the damping coefficient, and $G[n]$ is the energy density functional (which includes internal kinetic energy as well as exchange and correlation terms), which is a function of $n$, the free electron (number) density. We only account for the internal kinetic energy and as such, we use the Thomas-Fermi functional where the functional derivative is given as

$$
\frac{\delta G}{\delta n} =\frac{\hbar^2}{2m}(3\pi^2)^{2/3}n^{2/3}.
$$

For use in electrodynamics, we need to convert the velocity field to the current density $\textbf{J} = -en\textbf{v}$, where in general $n$ and $\textbf{v}$ are space and time-dependent.

With all of this, the full hydrodynamics formula is given by

$$
\frac{\partial \textbf{J}}{\partial t} =  \frac{\textbf{J}}{ne}\nabla\cdot \textbf{J}+\frac{1}{ne}(\textbf{J}\cdot\nabla)\textbf{J}-\frac{\textbf{J}}{en^2}(\textbf{J}\cdot\nabla)n+\frac{e^2n}{m}(\textbf{E}-\frac{\textbf{J}}{ne}\times\textbf{B})-\gamma\textbf{J}+\frac{e\hbar^2}{3m^2}(3\pi^2 n)^{2/3}\nabla n
$$

where

$$
n = -\frac{1}{e}\nabla\cdot \textbf{D}.
$$

These two formulas are sufficient for solving the hydrodynamic model.

## Numerical solution

The issue with this model is that it cannot be solved explicitly using Finite Differences. The following expression

$$
\frac{\partial \textbf{J}}{\partial t}\approx
\frac{\textbf{J}^{n+1/2}-\textbf{J}^{n-1/2}}{\Delta t}=\textbf{F}^n
$$



shows that we require every non-time-derivative term to be centered at $n$. Since we only know $J$ at $n-1/2$, we need to get creative. Inspired by [10.1109/JMMCT.2017.2751553], we solve the Hydrodynamic model using finite differences in two steps. Before we go to those two steps, let's start with finite differences for the time-derivative

 $$
\frac{\textbf{J}^{n+1/2}-\textbf{J}^{n-1/2}}{\Delta t} -\frac{\gamma}{2}(\textbf{J}^{n+1/2}+\textbf{J}^{n-1/2})\approx  \frac{e^2n}{m}\textbf{E}^n+\frac{\textbf{J}^n}{ne}\nabla\cdot \textbf{J}^n+\frac{1}{ne}(\textbf{J}^n\cdot\nabla)\textbf{J}^n-\frac{\textbf{J}^n}{en^2}(\textbf{J}^n\cdot\nabla)n-\frac{e^2n}{m}\frac{\textbf{J}^n}{ne}\times\textbf{B}^n+\frac{e\hbar^2}{3m^2}(3\pi^2 n)^{2/3}\nabla n
$$

which becomes

$$
\textbf{J}^{n+1/2} (\frac{1}{\Delta t}+\frac{\gamma}{2}) = \textbf{J}^{n-1/2} (\frac{1}{\Delta t}-\frac{\gamma}{2}) + \frac{e^2n}{m}\textbf{E}^n+\frac{\textbf{J}^n}{ne}\nabla\cdot \textbf{J}^n+\frac{1}{ne}(\textbf{J}^n\cdot\nabla)\textbf{J}^n-\frac{\textbf{J}^n}{en^2}(\textbf{J}^n\cdot\nabla)n-\frac{e^2n}{m}\frac{\textbf{J}^n}{ne}\times\textbf{B}^n+\frac{e\hbar^2}{3m^2}(3\pi^2 n)^{2/3}\nabla n
$$

The two steps become:

1.) We solve the (nonlocal) Drude model

$$
\textbf{J}_{Drude}^{n+1/2} (\frac{1}{\Delta t}+\frac{\gamma}{2}) = \textbf{J}^{n-1/2} (\frac{1}{\Delta t}-\frac{\gamma}{2}) + \frac{e^2n}{m}\textbf{E}^n+\frac{e\hbar^2}{3m^2}(3\pi^2 n)^{2/3}\nabla n
$$

which can be solved explicitly and used as an approximation for $\textbf{J}^{n+1/2}$. We then approximate $J^n \approx (J_{Drude}^{n+1/2} + J_{Full}^{n-1/2})/2$

2.) We then calculate a nonlinear correction:

$$
\textbf{J}_{corr}^{n+1/2} (\frac{1}{\Delta t}+\frac{\gamma}{2}) = \frac{\textbf{J}^n}{ne}\nabla\cdot \textbf{J}^n+\frac{1}{ne}(\textbf{J}^n\cdot\nabla)\textbf{J}^n-\frac{\textbf{J}^n}{en^2}(\textbf{J}^n\cdot\nabla)n-\frac{e^2n}{m}\frac{\textbf{J}^n}{ne}\times\textbf{B}^n.
$$

and the final current density is:
$$
\textbf{J}^{n+1/2}=\textbf{J}_{corr}^{n+1/2}+\textbf{J}_{Drude}^{n+1/2}.
$$

All finite differences are calculated using central differences. Linear interpolation between adjacent Yee cells is used when needed.

## How to use

Select a hydrodynamic material (codes 300 to 349: 300 is gold, 301 silver, 302 the material of
`[custom_material]`; see [Materials](Materials.md)) and the hydrodynamic electron model:

```ini
[geometry]
shape = sphere
material = 300

[electrons]
model = hydrodynamic

[hydrodynamic]
pressure = true          ; the quantum pressure (functional) term
convection = true        ; the convective term
magnetic = true          ; the magnetic Lorentz-force term
split_nonlinear = false  ; second-order perturbative model
```

| Setting | Model |
|---|---|
| all three terms on (default) | the full hydrodynamic model: best for small particles (diameter < 20 nm) with a small cell size (dx < 0.5 nm) |
| `pressure = false` | without the pressure term, $\nabla\frac{\delta G[n]}{\delta n} = 0$: the pressure only matters in small particles and at sharp features, and is otherwise not worth the computation |
| `convection`, `magnetic` | switch the two other nonlinear terms individually, to see what each contributes |
| `split_nonlinear = true` | the second-order perturbative model: the linear and the nonlinear fields are propagated separately, which makes second-order nonlinearities easy to analyse. Assumes $n = n_0$. Combine with `pressure = false` to drop the pressure term |

`[electrons] diffusion` and `fermi_velocity` scale the diffusion coefficient and the Fermi
velocity of the model.

For a second-harmonic spectrum analyse the pump and its harmonic:

```ini
[spectrum]
lambda_min = 800e-9
lambda_max = 800e-9
points = 1
harmonic_band = second
```

The peak field matters: set `[source] peak_field` [V/m] and use `signal = gaussian`.

A run with the hydrodynamic model cannot be continued from a state file yet
([Stopping and Restarting](../tutorial/Restart.md)).
