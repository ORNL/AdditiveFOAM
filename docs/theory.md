---
title: Theory
nav_order: 3
usemathjax: true
---

# Theory

---

## Governing equations

#### Mass and Momentum Conservation

$$
\nabla \cdot (\mathbf{u}) = 0
$$

$$
\rho\left(\frac{\partial \mathbf{u}}{\partial t} + \mathbf{u} \cdot \nabla \mathbf{u}\right) = \nabla \cdot (\mu \nabla \mathbf{u}) - \nabla p +  \rho _{k} \mathbf{g} - D\mathbf{u}
$$

where $\rho$ denotes density, $t$ denotes time, $\mathbf{u}$ denotes the mixture velocity vector, $\mu$ denotes the dynamic viscosity, and $p$ denotes pressusure. The source terms in the momentum equation represnt bouyancy induced from density variations due to temperature changes in the fluid and drag induced from coalesced solid.
 $\rho_{k}$ denotes the temperature-dependent density.

#### Energy Conservation
The heat equation without phase change is given by:

$$
\rho c_{p} \frac{\partial T}{\partial t} + \rho c_{p} \nabla \cdot \left( \mathbf{u} T\right) = \nabla \cdot \left(k\nabla T\right) + Q,
$$

where $$c_p$$ is the specific heat at constant pressure, $$T$$, is the
temperature, $$k$$ is the thermal conductivity, and $$Q$$ is the volumetric heat
source term.

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
146    // calculate the surface normal gradient on the patch
147    return
148    (
149        transform(I - sqr(nHat), pif) - pif
150     + coeff_*transform(I - sqr(nHat), tGrad) / this->patch().deltaCoeffs()
151   )*this->patch().deltaCoeffs();
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
