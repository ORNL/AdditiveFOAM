---
layout: versioned
title: Solver and Time Controls
parent: User Guide
nav_order: 8
permalink: /docs/solver-controls/
redirect_from:
  - /reference/solver-controls/
usemathjax: true
---

# Solver and Time Controls

This section controls four parts of the solution procedure: nonlinear temperature–phase correction, temporal discretization, optional PIMPLE fluid-flow coupling, and adaptive time-step selection. The controls are divided between the `PIMPLE` and linear-solver dictionaries in `system/fvSolution`, the `ddtSchemes` dictionary in `system/fvSchemes`, and the time controls in `system/controlDict`. A thermal-only case uses `nOuterCorrectors 0`; a coupled melt-pool-flow case uses a positive outer-corrector count.

## Thermodynamic coupling

These controls determine how many temperature–solid-fraction corrections are performed within each CFD time step. Keep `nThermoCorrectors` large enough for the reported maximum solid-fraction update to fall below `thermoTolerance`; the solver log reports that update for every correction. `Tmax` limits the temperature accepted by the implicit solve, while `explicitSolve` permits forward Euler only when the instantaneous diffusion criterion is below one.

AdditiveFOAM reads these entries from the `PIMPLE` dictionary in `system/fvSolution`:

| Entry | Default | Meaning |
|---|---:|---|
| `nThermoCorrectors` | `20` | Maximum temperature–solid-fraction iterations per time step |
| `thermoTolerance` | `1e-6` | Maximum solid-fraction update for convergence |
| `Tmax` | Unlimited | Optional implicit upper-temperature limiter |
| `explicitSolve` | `false` | Request explicit Euler when the diffusion criterion permits it |

```foam
PIMPLE
{
    momentumPredictor        no;
    nOuterCorrectors         0;
    nCorrectors              1;
    nNonOrthogonalCorrectors 0;
    pRefCell                 0;
    pRefValue                0;

    nThermoCorrectors        20;
    thermoTolerance          1e-8;
    Tmax                     4000;
    explicitSolve            false;
}
```

The thermodynamic loop performs at least two iterations and accepts convergence when

$$\max_c\left|f_{s,c}^{(m+1)}-f_{s,c}^{(m)}\right|<\varepsilon_f,$$

where $$c$$ indexes cells, $$m$$ is the thermodynamic-correction index, $$f_{s,c}^{(m)}$$ is the iterated solid fraction, and $$\varepsilon_f$$ is `thermoTolerance`; `nThermoCorrectors` bounds the iteration count.

## Time discretization

Select the temperature time scheme under `ddtSchemes` in `system/fvSchemes`. Use `Euler` for first-order implicit stepping, `backward` for second-order backward differencing, or `CrankNicolson` for an off-centred second-order scheme. Other schemes cause a fatal error. When the explicit path is active, AdditiveFOAM uses forward Euler regardless of the configured scheme. See [temperature and phase coupling]({{ '/docs/theory/#temperature-and-phase-coupling' | relative_url }}) for the coefficients used by each scheme.

## Flow coupling

AdditiveFOAM solves velocity and pressure with OpenFOAM's PIMPLE algorithm when `nOuterCorrectors` is positive. Each outer correction solves the momentum equation for `U`, applies the configured pressure corrections to `p_rgh`, and updates the velocity–pressure coupling.

{: .custom }
Set `nOuterCorrectors 0` for a thermal-only simulation. Set `nOuterCorrectors` to `1` or greater to enable fluid flow through the PIMPLE momentum–pressure loop. If there is no liquid in the domain, AdditiveFOAM skips that loop because the velocity is constrained to zero in fully solid material.

A thermal-only configuration is

```foam
PIMPLE
{
    momentumPredictor        no;
    nOuterCorrectors         0;
    nCorrectors              1;
    nNonOrthogonalCorrectors 0;
}
```

A coupled thermal-fluid configuration is, for example,

```foam
PIMPLE
{
    momentumPredictor        yes;
    nOuterCorrectors         1;
    nCorrectors              2;
    nNonOrthogonalCorrectors 0;
}
```

| Entry | Effect |
|---|---|
| `momentumPredictor` | Solves a momentum predictor before the pressure-correction loop when enabled |
| `nOuterCorrectors` | Number of complete momentum–pressure coupling loops; zero disables flow |
| `nCorrectors` | Pressure corrections performed within each outer loop |
| `nNonOrthogonalCorrectors` | Additional pressure corrections for mesh non-orthogonality |

### Linear solvers

`U` and `p_rgh` use the configured OpenFOAM linear solvers. `T` can be matched with a regular-expression solver entry:

```foam
solvers
{
    "T.*"
    {
        solver         PBiCGStab;
        preconditioner DILU;
        tolerance      1e-12;
        relTol         0;
    }
}
```

In this entry, `PBiCGStab` is the temperature-equation linear solver, `DILU` is its preconditioner, `tolerance` is the absolute residual stopping criterion, and `relTol 0` disables termination based on a relative residual reduction. Separate `U` and `p_rgh` entries are required when fluid flow is enabled.

## Adaptive time stepping

Set `adjustTimeStep yes` in `system/controlDict` to let AdditiveFOAM select the next time step from fluid advection, phase change, thermal diffusion, and scan-path event times. The solver prints the current Courant, thermo, and diffusion numbers and the accepted `deltaT`; compare those values with `maxCo`, `maxAlphaCo`, `maxDi`, and `maxDeltaT` to identify which limit controls the run.

```foam
adjustTimeStep yes;
maxCo          0.5;
maxDi          1;
maxAlphaCo     1;
maxDeltaT      1e-5;
```

For the current time step $$\Delta t^n$$, define the maximum phase-change number

$$A^n=\max_c\left(
\frac{|(\partial T/\partial t)_c|\Delta t^n}
{T_{\mathrm{liquidus}}-T_{\mathrm{solidus}}}
\right),$$

where $$n$$ is the current time-level index, $$T$$ is temperature, and $$T_{\mathrm{solidus}}$$ and $$T_{\mathrm{liquidus}}$$ bound the phase-change interval. The maximum finite-volume diffusion number is

$$D^n=\max_c\left[
\frac{\Delta t^n}{V_c\rho c_{p,c}}
\sum_{f\in\partial c}|\mathbf{S}_f|k_f\delta_f
\right].$$

In these definitions, $$A^n$$ is the maximum phase-change number, $$D^n$$ is the maximum diffusion number, $$c$$ indexes cells, and $$f\in\partial c$$ indexes the faces bounding cell $$c$$. The quantities $$V_c$$, $$c_{p,c}$$, and $$\rho$$ are cell volume, cell heat capacity, and density; $$\mathbf S_f$$ is the outward face-area vector, $$k_f$$ is face thermal conductivity, and $$\delta_f$$ is the inverse cell-centre spacing normal to the face. With fluid Courant number $$C^n$$, the unrestricted factor is

$$r=\min\left(
\frac{C_{\max}}{C^n+\epsilon},
\frac{A_{\max}}{A^n+\epsilon},
\frac{D_{\max}}{D^n+\epsilon}
\right).$$

Here $$C_{\max}$$, $$A_{\max}$$, and $$D_{\max}$$ are `maxCo`, `maxAlphaCo`, and `maxDi`; $$\epsilon$$ is a small positive denominator safeguard and $$r$$ is the unrestricted time-step ratio.

The applied growth factor and provisional step are

$$g=\min[r,1+0.1r,1.2],
\qquad
\Delta t^{n+1}=\min(g\Delta t^n,\Delta t_{\max}).$$

Here $$g$$ is the applied growth factor and $$\Delta t_{\max}$$ is `maxDeltaT`.

When `hitPathIntervals true`, the source model reduces this provisional value so that a step endpoint coincides with the next scan-path event. With `explicitSolve true`, forward Euler is selected only when the instantaneous $$D^n<1$$; otherwise the configured implicit scheme is used.
