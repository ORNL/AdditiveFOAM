---
title: Theory
nav_order: 3
usemathjax: true
---

# Theory

---

## Governing equations

The transport phenomena in `AdditiveFOAM` is modeled using continuum mixture relationships for mass, momentum, and energy.

#### Mass and Momentum Conservation
The conservation equations for mass and momentum are modeled under the assumption of incompressibility:

$$
\nabla \cdot (\mathbf{u}) = 0
$$

$$
\rho\left(\frac{\partial \mathbf{u}}{\partial t} + \mathbf{u} \cdot \nabla \mathbf{u}\right) = \nabla \cdot (\mu \nabla \mathbf{u}) - \nabla p +  \rho_k \mathbf{g} - D\mathbf{u}
$$

where $$\rho$$ denotes density, $$t$$ denotes time, $$\mathbf{u}$$ denotes the mixture velocity vector, $$\mu$$ denotes the dynamic viscosity, and $$p$$ denotes pressure. The source terms in the momentum equation represent buoyancy induced from density variations due to temperature changes in the liquid and drag induced from coalesced solid. The temperature-dependent density is calculated as:

$$
\rho_k = \rho \left[1 - \beta \left(T - T_{\text{liq}}\right)\right]
$$

where $$\beta$$ is the coefficient of thermal expansion, $$T$$ is temperature, and $$T_{\text{liq}}$$ is the alloy liquidus temperature. The drag induced from the coalesced solid is calculated from the Kozeny-Carman relationship:

$$
D = \mu \frac{180}{\lambda^{2}}  \frac{f_s^{2}}{(1 - f_s)^{3}}
$$

where $$\lambda$$ is the characteristic length scale for mushy zone drag and $$f_s$$ is the solid volme fraction.

{: .custom }
The value for $$\lambda$$ is chosen for numerical convenience to damp the flow in the rigid mush. We recommend using a value of 10e-6 microns

#### Energy Conservation
The conservation equation for energy is modeled in terms of temperature:

$$
\rho c_{p} \frac{\partial T}{\partial t} + \rho c_{p} \nabla \cdot \left( \mathbf{u} T\right) = \nabla \cdot \left(k\nabla T\right) + \rho L_{f} \frac{\partial f_s}{\partial t} + Q,
$$

where $$c_p$$ is the specific heat at constant pressure, $$k$$ is the thermal conductivity, $$L_f$$ is the latent heat of fusion, and $$Q$$ is the volumetric heat source term.

## Boundary Conditions
In addition to the standard boundary conditions provided in OpenFOAM, we provide two custom boundary conditions necessary for simulating AM processes.
The ```maragoni``` boundary condition is provided to calculate the shear stress induced by surface tensions gradients at along a boundary. The condition equates the viscous stress to the Marangoni stress:
$$
\mu \frac{\partial \mathbf{u}}{\partial \mathbf{n}} = \frac{\partial \gamma}{\partial T}\frac{\partial T}{\partial \mathbf{\tau}}
$$

Its implementation is derived from the ```partialSlip``` boundary conditions available in OpenFOAM-10, which uses two basic functions to update the boundary condition. The first function is ```Foam::marangoniFvPatchVectorField::snGrad()``` which sets the surface-normal gradient of the velocity field at the specified boundary faces. In OpenFOAM, ```snGrad()``` is notably used in the discretization of the laplacian operator in the momentum equation. Mathematically, the return of the function call is

$$
[\mathbf{I} - \mathbf{n}\cdot\mathbf{n}] \cdot\mathbf{u}_{i} - \mathbf{u}_{i} + c[\mathbf{I} - \mathbf{n}\cdot\mathbf{n}] \cdot \nabla T_{i}
$$

and is represented in the source code as:
```cpp
// calculate the surface normal gradient on the patch
return
(
    transform(I - sqr(nHat), pif) - pif
  + coeff_*transform(I - sqr(nHat), tGrad) / this->patch().deltaCoeffs()
)*this->patch().deltaCoeffs();
```
The first two terms uses vector transormations to subtract out the component of the velocoity field normal to the surface which may be numerically nonzero, and the last term equations the surface-normal gradient to the product of the maragoni coefficient $c=\frac{1}{\mu}\frac{\partial \gamma}{\partial T}$ and the component of the temperature gradient that is tangential to the boundary.

The second function equated in the boundary condition is ```Foam::marangoniFvPatchVectorField::evaluate()```

$$
\mathbf{u}_{b} = [\mathbf{I} - \mathbf{n}\cdot\mathbf{n}] \cdot\mathbf{u}_{i}
$$

```
vectorField::operator=(transform(I - sqr(nHat), pif));
```

The ```mixedTemperature``` boundary condition is provided to calculate the mixed convective and radiative heat transport across a boundary. 
