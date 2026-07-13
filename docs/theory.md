---
layout: versioned
title: Governing Equations
parent: User Guide
nav_order: 3
permalink: /docs/theory/
usemathjax: true
---

# Governing Equations

This section defines the mathematical model solved by AdditiveFOAM: incompressible momentum conservation with Boussinesq buoyancy and mushy-zone drag, temperature transport with latent heat, the semi-implicit correction that couples temperature to a tabulated solid-fraction path, supported time discretizations, and phase-dependent thermal properties with irreversible powder consolidation. Configuration entries and numerical algorithms that act on these equations are documented in the subsequent User Guide sections.

## Continuum model

AdditiveFOAM solves incompressible mass and momentum conservation for the mixture velocity $$\mathbf{u}$$:

$$\nabla\cdot\mathbf{u}=0,$$

$$\rho\left(\frac{\partial\mathbf{u}}{\partial t}+\mathbf{u}\cdot\nabla\mathbf{u}\right)
=\nabla\cdot(\mu\nabla\mathbf{u})-\nabla p+\rho_k\mathbf{g}-D\mathbf{u}.$$

The Boussinesq density factor is

$$\rho_k=1-\beta(T-T_{\mathrm{liq}}),$$

and the mushy-zone resistance is based on a Kozeny–Carman form,

$$D=\mu\frac{180}{\lambda^2}\frac{f_s^2}{(1-f_s)^3}.$$

Here $$\mathbf{u}$$ is mixture velocity, $$t$$ is time, $$\rho$$ is reference density, $$\mu$$ is dynamic viscosity, $$p$$ is pressure, $$\mathbf{g}$$ is gravitational acceleration, and $$\nabla$$ is the spatial gradient operator. The temperature is $$T$$, $$T_{\mathrm{liq}}$$ is the liquidus, $$\beta$$ is the thermal-expansion coefficient, and $$\rho_k$$ is the Boussinesq buoyancy factor. The solid fraction is $$f_s$$, $$D$$ is the mushy-zone drag coefficient, and $$\lambda$$ is the dendrite-arm spacing.

## Temperature and phase coupling

The energy equation is written in temperature form:

$$\rho c_p\frac{\partial T}{\partial t}
+\rho c_p\nabla\cdot(\mathbf{u}T)
=\nabla\cdot(k\nabla T)+\rho L_f\frac{\partial f_s}{\partial t}+\dot{q}.$$

Here $$c_p$$ is specific heat capacity, $$k$$ is thermal conductivity, $$L_f$$ is latent heat of fusion, and $$\dot q$$ is volumetric heat input. The `thermoPath` table defines the solid-fraction function $$f_s=F(T)$$ and can represent Gulliver–Scheil, lever-rule, or linear behavior. AdditiveFOAM interpolates piecewise linearly between the supplied points and applies $$\operatorname{clip}_{[0,1]}(x)=\min[\max(x,0),1]$$.

### Semi-implicit correction from the tabulated path

Let $$m$$ denote the thermodynamic-correction index at the new time level. For a temperature inside the solidification interval, the solver brackets $$T^{(m)}$$ by lower and upper table points $$(T_l,F_l)$$ and $$(T_h,F_h)$$ and evaluates the local path slope

$$g^{(m)}=\left.\frac{\mathrm{d}F}{\mathrm{d}T}\right|_m
=\frac{F_h-F_l}{T_h-T_l}.$$

Rather than anchoring the linearization at a temperature that may not agree with the current iterated solid fraction, it maps $$f_s^{(m)}$$ back onto that segment:

$$T_0^{(m)}=T_l+\frac{f_s^{(m)}-F_l}{g^{(m)}}.$$

Thus the current phase state lies exactly on the local linear model. The next solid-fraction estimate is

$$f_s^{(m+1)}=
\operatorname{clip}_{[0,1]}
\left[f_s^{(m)}+g^{(m)}\left(T^{(m+1)}-T_0^{(m)}\right)\right].$$

For a linear multistep time derivative, only the new-time value changes during the nonlinear correction. If $$a_0$$ is the new-time derivative coefficient and $$\Delta t$$ is the current time-step duration, then

$$\left.\frac{\partial f_s}{\partial t}\right|^{(m+1)}
\approx
\left.\frac{\partial f_s}{\partial t}\right|^{(m)}
+\frac{a_0}{\Delta t}g^{(m)}
\left(T^{(m+1)}-T_0^{(m)}\right).$$

Substitution into the energy equation gives the temperature solve

$$\mathcal{L}_T\!\left(T^{(m+1)}\right)
=\rho L_f\left[
\left.\frac{\partial f_s}{\partial t}\right|^{(m)}
+\frac{a_0}{\Delta t}g^{(m)}
\left(T^{(m+1)}-T_0^{(m)}\right)
\right]+\dot q,$$

where $$\mathcal{L}_T$$ contains sensible-energy storage, advection, and diffusion. Moving the term proportional to $$T^{(m+1)}$$ to the matrix gives

$$\mathcal{L}_T\!\left(T^{(m+1)}\right)
-\rho L_f\frac{a_0}{\Delta t}g^{(m)}T^{(m+1)}
=\rho L_f\left[
\left.\frac{\partial f_s}{\partial t}\right|^{(m)}
-\frac{a_0}{\Delta t}g^{(m)}T_0^{(m)}
\right]+\dot q.$$

Within the usual decreasing solid-fraction path, $$g^{(m)}\leq0$$, so the latent-heat Jacobian adds a non-negative contribution to the temperature-equation diagonal. This is the semi-implicit stabilization: the tabulated path remains explicit data, while its local slope enters the linear system implicitly.

The new temperature is used in the clipped phase update above. Iteration stops after at least two corrections when

$$\max_{\text{cells}}\left|f_s^{(m+1)}-f_s^{(m)}\right|
<\texttt{thermoTolerance},$$

or when `nThermoCorrectors` is reached. Outside the tabulated solidification interval, limiting branches drive inconsistent phase states back toward the liquidus or solidus.

### Time-integration methods

The new-time coefficient used in the correction is

| Scheme | $$a_0$$ | Behavior |
|---|---:|---|
| `Euler` | $$1$$ | First-order implicit Euler |
| `backward` | $$1+\dfrac{\Delta t}{\Delta t+\Delta t_0}$$ | Second-order backward differencing with variable time steps; $$a_0=3/2$$ for equal steps |
| `CrankNicolson` | $$1+\psi$$ | Crank–Nicolson with off-centering coefficient $$\psi$$; $$\psi=1$$ is fully centred |

Here $$\Delta t_0$$ is the preceding time-step duration and $$\psi$$ is the OpenFOAM Crank–Nicolson off-centering coefficient.

Select the scheme in `system/fvSchemes`:

```foam
ddtSchemes
{
    default backward;
}
```

and configure the nonlinear correction in `system/fvSolution`:

```foam
PIMPLE
{
    nThermoCorrectors 20;
    thermoTolerance   1e-8;
    explicitSolve     false;
}
```

`explicitSolve true` requests explicit sensible-energy transport when its diffusion-number check passes and forces forward Euler regardless of `fvSchemes`. Use `false` to guarantee implicit `Euler`, `backward`, or `CrankNicolson` transport.

### Phase properties and powder transition

Let $$f_p$$ be the `alpha.powder` field and $$f_s$$ the solid fraction obtained from `thermoPath`. The liquid fraction is

$$f_l=1-f_s.$$

For either thermal conductivity or heat capacity, denoted generically by $$\psi\in\{k,c_p\}$$, AdditiveFOAM evaluates the mixture property as

$$\psi_{\mathrm{eff}}
=(1-f_p)\left[f_s\psi_s(T_s)+f_l\psi_l(T_l)\right]
+f_p\psi_p(T_s),$$

where $$\psi_s$$, $$\psi_l$$, and $$\psi_p$$ are the configured solid, liquid, and powder property polynomials. Their evaluation temperatures are limited to

$$T_s=\operatorname{clip}(T,300\ \mathrm{K},T_{\mathrm{solidus}}),$$

$$T_l=\operatorname{clip}(T,T_{\mathrm{liquidus}},2T_{\mathrm{liquidus}}).$$

Here $$f_p$$, $$f_s$$, and $$f_l$$ are the powder, solid, and liquid fractions; $$T_s$$ and $$T_l$$ are the temperatures used to evaluate the solid/powder and liquid polynomials, respectively. The subscript on $$\psi_i$$ identifies phase $$i\in\{s,l,p\}$$.

Powder consolidation is represented by the irreversible update

$$f_p^{\,n+1}=
\begin{cases}
0, & f_s^{\,n+1}<1,\\
f_p^{\,n}, & f_s^{\,n+1}=1.
\end{cases}$$

Here $$n$$ and $$n+1$$ denote the old and new time levels.

Thus a cell retains its initialized powder fraction while it remains completely solid. As soon as the cell begins melting, its powder fraction is set to zero and the cell thereafter uses the consolidated solid–liquid mixture, even if it later resolidifies.
