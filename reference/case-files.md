---
layout: versioned
title: Case Files and Materials
parent: User Guide
nav_order: 8
permalink: /docs/case-files/
redirect_from:
  - /reference/case-files/
usemathjax: true
---

# Case Files and Materials

This section describes the files needed to construct an AdditiveFOAM case: initial fields, bundled and user-defined materials, temperature-dependent phase properties, the tabulated temperature–solid-fraction relation, scan-path input, and the OpenFOAM dictionaries that control the mesh and solution. AdditiveFOAM currently provides ready-to-use configurations for IN625, SS316L, and AlSi10Mg under `$ADDITIVEFOAM_ETC/materials`.

## Case structure

An AdditiveFOAM case follows the OpenFOAM directory structure:

```text
case/
├── 0/
│   ├── T
│   ├── U
│   ├── p_rgh
│   └── alpha.powder
├── constant/
│   ├── transportProperties
│   ├── heatSourceDict
│   ├── scanPath
│   └── dynamicMeshDict
└── system/
    ├── blockMeshDict
    ├── controlDict
    ├── fvSchemes
    ├── fvSolution
    └── decomposeParDict
```

The `0/` directory defines initial and boundary values. `constant/transportProperties` defines the material and thermodynamic path; `constant/heatSourceDict` selects the moving sources, absorption models, source integration, and refinement model; each scan-path file defines source motion and power. `constant/dynamicMeshDict` configures OpenFOAM refinement and load balancing when those features are used.

Under `system/`, `blockMeshDict` defines the initial mesh, `controlDict` defines run time and Function Objects, `fvSchemes` selects discretizations, `fvSolution` selects linear solvers and PIMPLE controls, and `decomposeParDict` defines the parallel decomposition and redistribution method.

## Initial fields

| Field | Read behavior | Notes |
|---|---|---|
| `T` | Required | Temperature with initial and boundary values |
| `U` | Optional | Defaults to zero; required for meaningful fluid-flow initialization |
| `p_rgh` | Optional | Defaults to zero |
| `alpha.powder` | Optional | Defaults to zero; powder is irreversibly removed once melting begins |
| `alpha.solid` | Not read | Calculated from `T` and `thermoPath` |

`T` must define every mesh patch. `U` and `p_rgh` are required for a fluid-flow simulation but default to zero when absent. `alpha.powder` identifies the initially unconsolidated region; it defaults to zero for a fully consolidated domain. The solver creates `alpha.solid` from the initial temperature and writes it as a calculated field.

## Materials

Sourcing `etc/bashrc` sets `ADDITIVEFOAM_ETC`. AdditiveFOAM currently supplies these material configurations under `$ADDITIVEFOAM_ETC/materials/`:

| Configuration | Material | Solidus (K) | Liquidus (K) | Tutorials |
|---|---|---:|---:|---|
| `IN625.cfg` | Inconel 625 | 1410 | 1620 | AMB2018-02-B, multi-beam, multi-layer PBF |
| `SS316L.cfg` | 316L stainless steel | 1471 | 1709 | nLight AFX |
| `AlSi10Mg.cfg` | AlSi10Mg aluminum alloy | 850 | 870 | tabulated source |

Include one configuration from `constant/transportProperties`:

```foam
#include "$ADDITIVEFOAM_ETC/materials/IN625.cfg"
```

Each material file supplies:

- `solid`, `liquid`, and `powder` polynomial coefficients for `kappa` and `Cp`;
- `rho`, `mu`, `beta`, `DAS`, and `Lf`;
- `emissivity` for `mixedTemperature` and `dSigmadT` for `marangoni`;
- the material `thermoPath`.

A material not listed above is defined by providing the same entries directly in `transportProperties` or in another included dictionary.

## `transportProperties`

Use a bundled material include when one of the supplied alloys applies. For another material, define the same phase polynomials and constant properties directly in `constant/transportProperties`; these values control heat conduction, sensible-energy storage, buoyancy, viscous flow, mushy-zone resistance, and latent heat.

Each phase supplies quadratic polynomial coefficients in ascending powers of temperature:

```foam
solid
{
    kappa (a0 a1 a2);
    Cp    (b0 b1 b2);
}

liquid
{
    kappa (a0 a1 a2);
    Cp    (b0 b1 b2);
}

powder
{
    kappa (a0 a1 a2);
    Cp    (b0 b1 b2);
}
```

For each phase $$i\in\{s,l,p\}$$, the corresponding coefficient list defines

$$\psi_i(T)=c_{i,0}+c_{i,1}T+c_{i,2}T^2,
\qquad \psi\in\{k,c_p\}.$$

Here $$i=s,l,p$$ denotes the solid, liquid, or powder phase; $$\psi_i$$ is either thermal conductivity $$k_i$$ or heat capacity $$c_{p,i}$$; $$T$$ is temperature; and $$(c_{i,0},c_{i,1},c_{i,2})$$ denotes the applicable `kappa` or `Cp` coefficient list.

The solid and powder polynomials are evaluated at $$T_s=\operatorname{clip}_{[300\ \mathrm{K},\,T_{\mathrm{solidus}}]}\!\left[T\right]$$, and the liquid polynomial at $$T_l=\operatorname{clip}_{[T_{\mathrm{liquidus}},\,2T_{\mathrm{liquidus}}]}\!\left[T\right]$$. The [phase-properties section]({{ '/docs/theory/#phase-properties-and-powder-transition' | relative_url }}) defines their mixture weights and the irreversible powder transition.

The phase conductivities and heat capacities vary with temperature through the polynomials above. The remaining transport properties are constants within a simulation:

| Entry | SI units | Role in the model |
|---|---|---|
| `rho` | kg/m³ | Reference density in momentum, sensible energy, and latent energy |
| `mu` | Pa·s | Liquid dynamic viscosity and mushy-zone drag coefficient |
| `beta` | 1/K | Boussinesq thermal-expansion coefficient |
| `DAS` | m | Dendrite-arm spacing or equivalent mushy-zone drag length $$\lambda$$ |
| `Lf` | J/kg | Latent heat of fusion $$L_f$$ |

The shared entries `dSigmadT` and `emissivity` provide default material values for `marangoni` and `mixedTemperature`. A patch-local value overrides the shared value.

## `thermoPath`

Set `thermoPath` to the temperature–solid-fraction relation for the selected alloy. The solver uses it to initialize `alpha.solid`, update phase fraction during every thermodynamic correction, evaluate latent heat, identify the solidus and liquidus, and provide default isovalues to heat-source and solidification outputs.

The preferred 2.0 form is an entry inside `constant/transportProperties`:

```foam
thermoPath
(
    (1410 1)
    (1620 0)
);
```

For compatibility, if that entry is absent, AdditiveFOAM reads the same list from the legacy `constant/thermoPath` file. At least two points are required, and the path must contain exact entries with solid fraction `1` and `0` so the solidus and liquidus can be inferred. For adjacent table points $$(T_i,f_i)$$ and $$(T_{i+1},f_{i+1})$$,

$$f_s(T)=\operatorname{clip}_{[0,1]}\!\left[
f_i+\frac{f_{i+1}-f_i}{T_{i+1}-T_i}(T-T_i)
\right].$$

Here $$i$$ is the table-interval index, $$T_i$$ and $$T_{i+1}$$ are consecutive temperatures, and $$f_i$$ and $$f_{i+1}$$ are their tabulated solid fractions.

The table may represent any tabulated temperature–solid-fraction relation, including Gulliver–Scheil, lever-rule, or linear behavior. Its liquidus is the default isovalue for transient source depth, ExaCA, and solidification data.

## Scan-path file

The first line is a header and is ignored. Each subsequent non-empty line contains:

```text
mode x y z power parameter
```

| Mode | Position | Parameter |
|---|---|---|
| `0` | Segment endpoint in metres | Translation speed in m/s |
| `1` | Fixed spot position in metres | Dwell time in seconds |

Power is in watts. For a mode-0 row connecting endpoints $$\mathbf{x}_i$$ and $$\mathbf{x}_{i+1}$$ at speed $$v_i$$, the duration and position are

$$\Delta t_i=\frac{\|\mathbf{x}_{i+1}-\mathbf{x}_i\|}{v_i},
\qquad
\mathbf{x}(t)=\mathbf{x}_i+
\frac{t-t_i}{\Delta t_i}(\mathbf{x}_{i+1}-\mathbf{x}_i).$$

Here $$i$$ is the scan-path interval index, $$\mathbf{x}_i$$ and $$\mathbf{x}_{i+1}$$ are its endpoints, $$v_i$$ is its prescribed speed, $$t_i$$ is its start time, and $$\Delta t_i$$ is its duration.

For mode 1, $$\mathbf{x}(t)=\mathbf{x}_i$$ and the parameter is $$\Delta t_i$$. A zero-power interval contributes no heat, and zero-duration intervals are omitted. Path time begins at zero; the interval search is bidirectional so restart times earlier or later than the current cached interval are supported.

`hitPathIntervals` defaults to `true`, causing adaptive time stepping to land on path event boundaries. The solver stops evaluating a beam after the last positive-power segment; trailing zero-power motion does not extend active heating.

## Time-step controls

The [solver and time-controls section]({{ '/docs/solver-controls/#adaptive-time-stepping' | relative_url }}) defines the Courant, phase-change, and diffusion criteria and the resulting update equation for $$\Delta t$$.
