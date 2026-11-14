# Numerical Dispersion 

Due to the finitely spaced grid the simulated electromagnetic wave is subject to numerical dispersion. 
In the Total-Field / Scattered - Field approach this dispersion has to be taken into account in order to minimize the 
leakage of power into the scattered field region from the incident wave. 
For example, simulating an empty domain the wave should not be able to "escape" the TFSF box.
If we do not correct for numerical dispersion this might happen, as we assume a slightly wrong wave vector and thus a wrong phase velocity of the wave.

In Taflove, the numerical dispersion relation is derived for a plane wave. However, we are always using an envelope function, like e.g. the Gaussian. 
Therefore the dispersion relation is not 100% correct. We would have to derive the correct relation for our various envelope functions.
However, we just implemented the plane wave dispersion relation so far. 

You can probably do this implementation also for signals with envelopes, but you might have to update the numerical wave vector at every time step.

# Implementation 
First of all, we only have to correct for dispersion when the implementation of oblique incidence is used. For the standard normal incidence we have a 1D simulation that is run before the actual 3D simulation to account for the numerical dispersion.

The correction is controlled by `[wave_superposition] dispersion_correction` (default `true`).

If activated, the correction, as described below, is performed in `src/sources/plane_waves.cc`.

The numerical dispersion for a plane wave derived in Taflove is the following:

$$
\left[ \frac{1}{c \Delta t} sin\left(\frac{\omega \Delta t}{2}\right) \right]^2 = \left[\frac{1}{\Delta x} sin\left(\frac{\tilde{k}_x \Delta x }{2}\right)\right]^2 +  \left[\frac{1}{\Delta y} sin\left(\frac{\tilde{k}_y \Delta y }{2}\right)\right]^2 +  \left[\frac{1}{\Delta z} sin\left(\frac{\tilde{k}_z \Delta z }{2}\right)\right]^2
$$

We want to find the numerical wave vector $\tilde{k}$ in order to use it as the wave vector in the excitation functions (`src/sources/oblique.cc`). From the input file we know the direction of the wave vector $\hat{k} = (\hat{k}_x, \hat{k}_y, \hat{k}_z)$. Hence we can express the numerical wave vector as $ \tilde{k} = | \tilde{k}| \cdot \hat{k}$. We transformed the function into a relation that only depends on one variable, namely $|\tilde{k}|$. 

To solve for $|\tilde{k}|$ we can use the [`Newton-Raphson-Method`](https://en.wikipedia.org/wiki/Newton%27s_method). Therfore we define a new function $f(r)$, which is just the dispersion relation from above brought to one side.

$$
f(r) = - \left[ \frac{1}{c \Delta t} sin\left(\frac{\omega \Delta t}{2}\right) \right]^2 + \left[\frac{1}{\Delta x} sin\left(\frac{r \cdot \hat{k}_x \Delta x }{2}\right)\right]^2 +  \left[\frac{1}{\Delta y} sin\left(\frac{r \cdot \hat{k}_y \Delta y }{2}\right)\right]^2 +  \left[\frac{1}{\Delta z} sin\left(\frac{r \cdot \hat{k}_z \Delta z }{2}\right)\right]^2
$$

$|\tilde{k}|$ satisfies the numerical dispersion relation when $f(|\tilde{k}|) = 0$. For applying the Newton-Raphson method, we also need the derivation of $f(r)$.

$$
f'(r) = 
$$

We implemented the Newton-Raphson method in `src/sources/plane_waves.cc`. Currently, it is using a fixed number of iterations and a fixed numerical precision.  We could make these arguments later on.

The method returns a new absolute value of the wave vector $|\tilde{k}|$. THe new wavevector implies an updated phase velocity of the wave (smaller than $c_0$). 

$$
\tilde{c} = \frac{\omega}{|\tilde{k}|} = D \cdot c_0.
$$

We introduce here the dispersion factor $D$, which we can multiply by the speed of light to obtain th phase velocity corrected for dispersion.

$$
D = \frac{\tilde{c}}{c_0} = \frac{2 \pi f}{c_0 |\tilde{k}|}
$$

The dispersion factor is stored with the superposition source (`dispersion_factor`) and used in the excitation functions (`src/sources/oblique.cc`) to define the phase velocity.