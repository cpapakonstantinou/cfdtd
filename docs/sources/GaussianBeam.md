
# Gaussian Beam Theory and Code Usage

## Overview

In the code a Gaussian Beam can be simulated as a source for an FDTD simulation. This code originally derived from the paper 'Computation of tightly-focused laser beams in the FDTD method' [doi link here](https://doi.org/10.1364/OE.21.000087). The following, highlights the connections between the code and the paper. Note, some liberties were taken with the variables but the integrity of the equations and code remains. 
The last section shows how to use the source.

## Theory to Code
### The Theory Equations

Beginning with equation 8 in the [paper](https://doi.org/10.1364/OE.21.000087).

$$E(r)=\frac{n_m}{2\pi c} \sum{\alpha_n \hat{a} (\theta_n,\phi_n) \dot{E} ^\infty(\theta_n,\phi_n,t-\frac{n_m}{c} \hat{s} _n \cdot r)}$$

Replacing the field strength ($$E^\infty$$) with the definition given in equation 6.

$$E(r)=\frac{n_m}{2\pi c}\sum{\alpha_n \hat{a}(\theta_n,\phi_n)\frac{d}{dt} \sqrt{\frac{n_{\tilde{m}}}{n_m}}f\sqrt{\cos{\theta}}\psi(t-\frac{n_m}{c}\hat{s} _n\cdot r)\gamma_{gh}(x,y)}$$

Using the definition of gamma from equation 1 we acheive a mostly expanded form (expanded in all the ways that matter for connecting the pins).

$$E(r)=\frac{n_m}{2\pi c}\sum{\alpha_n \hat{a}(\theta_n,\phi_n)\frac{d}{dt} \sqrt{\frac{n_{\tilde{m}}}{n_m}}f\sqrt{\cos{\theta}}\psi(t-\frac{n_m}{c}\hat{s} _n\cdot r) H_g(\frac{\sqrt{2}x}{\omega_0})H_h(\frac{\sqrt{2}y}{\omega_0})e^{-\frac{x^2+y^2}{\omega_0^2}}}$$
### The Code Equations

In the plane-wave excitation of each face (`src/sources/oblique.cc`) the following are defined, here for E_z.

`To_z = (K * DZ - FOCUS_Z) * s_z;`
`Tm_z = ((K + 0.5) * DZ - FOCUS_Z) * s_z;`
`delay_z = (To_x + To_y + Tm_z) * N_RIFR[first_medium] / C0;`
`Ez_excit = u_z * signal_derivative(t_instant - d_z);`

Where u_z is the z component of the electric amplitude.
In mathematical notation these translate to the following respectively. Note 'dSig' is a function.

$$To_z=(I\cdot dz-f_z)s_z$$
$$Tm_z=((I+0.5)dz - f_z)s_z$$
$$delz=(To_x+To_y+Tm_z)\frac{n_m}{c_0}$$
$$E_z=u_z \cdot dSig(t-delz)$$

Combining these equations:

$$E_z=u_zdSig(t-[(I\cdot dz-f_z)s_z+To_y+((I+0.5)dz - f_z)s_z]\frac{n_m}{c_0}).$$

In `build_oblique_plane_wave` (`src/sources/gaussian_beam.cc`) $$u_z$$ is defined.

`u_z = AMPL * (cos_phi * cos_phi + cos_theta * sin_phi * sin_phi)`

Utilizing the default amplitude and exponent the code declares the following variables.
`exponent = FOCAL_LENGTH * sin_theta / BEAM_W0;`
`AMPL = AMPL_CONST * sqrt(cos_theta) * sin_theta * exp(-exponent * exponent);`

Where AMPL_CONST is defined as $$\frac{nf\Delta\theta\Delta\phi}{2\pi c_0}$$ in the simulation setup file. Translating these coded equations into mathematical notation and combining them we get the following.

$$u_z=\frac{nf}{2\pi c_0}\Delta\theta\Delta\phi\sqrt{\cos{\theta}}\sin{\theta}e^{-(\frac{f\sin{\theta}}{w_0})^2}(\cos^2{\phi}+\cos{\theta}\sin^2{\theta})$$

Plugging this into the electric component equation we achieve the fully expanded form.

$$E_z=\frac{nf}{2\pi c_0}\Delta\theta\Delta\phi\sqrt{\cos{\theta}}\sin{\theta}e^{-(\frac{f\sin{\theta}}{w_0})^2}(\cos^2{\phi}+\cos{\theta}\sin^2{\phi}) Sig(t-[(I\cdot dz-f_z)s_z+To_y+((I+0.5)dz - f_z)s_z]\frac{n_m}{c_0})$$

### The Theory-Code Parallels

Before drawing comparisons, first pointing out the differences in the equations that have thus far been derived. 1) The paper equation is a function of r while the code equation is only the z component of the electric field. 2) The paper equation sums over all possible values of $$\phi$$ and $$\theta$$ while the code equation comes from inside a double for loop over the angles $$\phi$$ and $$\theta$$. Additionally, the paper is focused on Hermite-Gaussian beams, but the code only focuses on Gaussian beams (the (0,0) mode of Hermite-Gaussian beams), therefore we only consider $$H_0=1$$, this change has been made for the following comparisons. The equations shown below are the same as derived above but grouped to highlight the similarities.

Paper:
$$E(r)=[{\frac{n_m f}{2\pi c}}]\sum{e^{-\frac{x^2+y^2}{\omega_0^2}}\sqrt{\cos{\theta}} \,\, [\alpha_n \sqrt{\frac{n_{\tilde{m}}}{n_m}} \hat{a}(\theta_n,\phi_n)]\frac{d}{dt} \psi(t-[\frac{n_m}{c}]\hat{s} _n\cdot r) }$$
Code:
$$E_z=[\frac{nf}{2\pi c_0}]e^{-(\frac{f\sin{\theta}}{w_0})^2}\sqrt{\cos{\theta}}[\Delta\theta\Delta\phi\sin{\theta}(\cos^2{\phi}+\cos{\theta}\sin^2{\phi})]dSig(t-((I\cdot dz-f_z)s_z+To_y+((I+0.5)dz - f_z)s_z)[\frac{n}{c_0}])$$

Listing the similarities moving left to right through the equations (left side of mathematical equality from paper, right side from code):
1. These constants are equivalent.
$$\frac{n_m f}{2\pi c} \equiv \frac{nf}{2\pi c_0}$$
2. The exponential uses different coordinate systems but is equivalent.
$$e^{-\frac{x^2+y^2}{\omega_0^2}}\equiv e^{-(\frac{f\sin{\theta}}{\omega_0})^2}$$
3. The square root of cosine theta is equivalent.
$$\sqrt{\cos{\theta}} \equiv \sqrt{\cos{\theta}}$$ 
4. The weighting and directional polarization align with the trigonometric quantities. 
$$\alpha_n \sqrt{\frac{n_{\tilde{m}}}{n_m}} \hat{a}_n (\theta_n,\phi_n) \equiv \Delta\theta\Delta\phi\sin{\theta}(\cos^2{\phi}+\cos{\theta}\sin^2{\phi})$$
5. And, the time derivative of the time shifted wavefront corresponds to the coded envelope produced by the time shifted Signal.
$$\frac{d}{dt} \psi(t-\frac{n_m}{c}\hat{s} _n \cdot r) \equiv dSig(t-\frac{n}{c_0}((I(dz)-f_z)s_z+To_y+((I+0.5)dz-f_z)s_z)$$

## Using the Code

`[source] type = gaussian_beam` builds the beam as a superposition of plane waves on a grid of
polar and azimuthal angles and injects them on the TF/SF box.

| `[gaussian_beam]` key | Meaning |
|---|---|
| `theta_steps`, `phi_steps` | number of intervals of the polar angle θ and of the azimuth φ; the beam is `theta_steps * phi_steps` plane waves |
| `theta_max` | [deg] largest polar angle (the numerical aperture) |
| `phi_max` | [deg] range of the azimuth; 360 for a round beam |
| `focal_factor` | focal length = `focal_factor` × the distance from the TF/SF corner to the structure origin |
| `waist_factor` | beam waist = `waist_factor` × sin(`theta_max`) × focal length |

- **Focal length:** a factor below 1 brings the figurative lens inside the domain (for an origin at
  the centre of the domain); above 1 moves it away from the domain.
- **Beam waist:** keep the waist in proportion to the focal length; an extreme `waist_factor` gives a
  warped intensity profile.
- **Angles:** `theta_max` slightly above 45° worked well. Much higher and the waves interfere
  constructively outside the beam; much lower and the focus stays weak. `phi_max` is measured from
  the x axis in the x-z plane.
- **Number of waves:** every wave is evaluated at every time step, so the run time grows with
  `theta_steps * phi_steps`. Too few waves give an asymmetric beam that does not converge; keep
  `phi_steps` a multiple of 4 for symmetry.
- The structure origin fixes the focus, so the beam cannot be moved without changing the domain.

### Example

```ini
[grid]
nx = 160
ny = 400
nz = 160
dx = 50e-9
dy = 50e-9
dz = 50e-9
time_steps = 2000

[parallel]
ranks_x = 2
ranks_y = 4
ranks_z = 2

[boxes]
tfsf_x = 10
tfsf_y = 10
tfsf_z = 10

[spectrum]
lambda_min = 0.5e-6
lambda_max = 1.5e-6
points = 31

[source]
type = gaussian_beam
signal = gaussian

[signal.gaussian]
center_wavelength = 1e-6
width = 25e-15

[gaussian_beam]
theta_steps = 15
phi_steps = 60
theta_max = 50
phi_max = 360
waist_factor = 0.22
focal_factor = 1.0

[output]
frames_xy = 15
frames_xz = 15
```

A hand-built list of waves does the same through `type = wave_superposition`
([Plane Wave Superposition](PlaneWaveSuperposition.md)).
