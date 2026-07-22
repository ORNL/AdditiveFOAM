---
layout: versioned
title: Function Objects
parent: User Guide
nav_order: 10
permalink: /docs/function-objects/
redirect_from:
  - /reference/function-objects/
usemathjax: true
---

# Function Objects

AdditiveFOAM provides three Function Objects for extracting melt-pool geometry and solidification conditions while a case runs:

| Function Object | Result |
|---|---|
| `meltPoolDimensions` | Time histories of isotherm length, width, and depth |
| `solidificationData` | Cell-centred cooling rate, temperature gradient, and isotherm velocity at downward temperature crossings |
| `ExaCA` | Melt time, solidification time, and cooling rate on a user-defined Cartesian sampling lattice for ExaCA |

Configure any of these under `functions` in `system/controlDict` and load `libadditiveFoamFunctionObjects.so`. `meltPoolDimensions` writes continuously during the simulation; `solidificationData` and `ExaCA` collect crossing events and write per-process data that is combined with the reconstruction utilities.

## Melt-pool dimensions

Use `meltPoolDimensions` when the quantities of interest are the overall length, width, and depth of one or more temperature isosurfaces. Add the object before running the case, choose `isoValues`, and use `scanPathAngle` to align reported length with the scan direction. The object writes one CSV time history per isovalue, which `plotDimensions` converts into a figure.

```foam
meltPoolDimensions
{
    libs          ("libadditiveFoamFunctionObjects.so");
    type          meltPoolDimensions;
    enabled       true;
    isoValues     (1410 1620);
    scanPathAngle 0;
}
```

`isoValues` defaults to the solidus and liquidus from `thermoPath`. For each internal or coupled face joining cell centres $$\mathbf{x}_o$$ and $$\mathbf{x}_n$$ whose temperatures bracket an isovalue $$T_I$$, the estimated isosurface point is

$$\mathbf{x}_I=\mathbf{x}_o+
\frac{T_I-T_o}{T_n-T_o}(\mathbf{x}_n-\mathbf{x}_o).$$

Here subscripts $$o$$ and $$n$$ denote the owner and neighbour cells, $$T_o$$ and $$T_n$$ are their temperatures, and $$\mathbf{x}_I$$ is the linearly interpolated position at which the temperature equals isovalue $$T_I$$.

For `scanPathAngle` $$\theta$$, each point is transformed as

$$x'=x\cos\theta+y\sin\theta,\qquad
y'=y\cos\theta-x\sin\theta,\qquad z'=z.$$

Here $$(x,y,z)$$ are the original coordinates, $$(x',y',z')$$ are the rotated coordinates, and $$\theta$$ is `scanPathAngle` converted from degrees to radians.

Physical-boundary face centres are included when either the boundary or adjacent-cell temperature is at least $$T_I$$. For the resulting global set $$\mathcal{S}_I$$, the reported dimensions are

$$L=\max_{\mathcal{S}_I}x'-\min_{\mathcal{S}_I}x',\quad
W=\max_{\mathcal{S}_I}y'-\min_{\mathcal{S}_I}y',\quad
D=\max_{\mathcal{S}_I}z'-\min_{\mathcal{S}_I}z'.$$

Here $$\mathcal S_I$$ is the global set of points associated with isovalue $$T_I$$, and $$L$$, $$W$$, and $$D$$ are its streamwise length, transverse width, and vertical depth.

The master rank writes one file per isovalue:

```text
postProcessing/meltPoolDimensions/<isovalue>.csv
```

Columns are `time(s),length(m),width(m),depth(m)`. Plot them with `plotDimensions`.

## Solidification data

Use `solidificationData` to record local solidification conditions at finite-volume cell centres. Define `box` around the region of interest and choose the crossing temperature with `isoValue`; the material liquidus is used by default. After a parallel run, combine the per-process event files with `reconstructSolidificationData`, then use `plotCET` to plot thermal gradient against isotherm velocity.

```foam
solidificationData
{
    libs    ("libadditiveFoamFunctionObjects.so");
    type    solidificationData;
    enabled true;
    box     (-1 -1 -1) (1 1 1);
    isoValue 1620;
}
```

`box` is required and `isoValue` $$T_I$$ defaults to the liquidus. A cell event is registered at time $$t^n$$ when

$$T_c^{n-1}>T_I,\qquad T_c^n\leq T_I.$$

At the cell centre, the output quantities are

$$R_c=\left|\frac{\partial T_c}{\partial t}\right|,\qquad
\mathbf{G}_c=(\nabla T)_c,\qquad
V_c=\frac{R_c}{\max(|\mathbf{G}_c|,\epsilon)}.$$

Here $$c$$ identifies the crossing cell, $$R_c$$ is cooling-rate magnitude, $$\mathbf G_c$$ is the cell-centred temperature gradient, $$V_c$$ is the estimated isotherm speed, and $$\epsilon$$ is a small positive safeguard for a zero gradient.

The CSV columns are `x,y,z,ts,cr,gx,gy,gz,v`; `ts` is $$t^n$$, `cr` is the stored magnitude $$R_c$$, and `v` is the local isotherm-speed estimate $$V_c$$. Each rank writes `solidificationData/data_<rank>.csv` at the end of the run. `reconstructSolidificationData` combines them into `solidificationData/solidification-data.csv`, and `plotCET` plots $$\lvert\mathbf{G}\rvert$$ against $$V$$.

## ExaCA thermal histories

Use `ExaCA` when exporting melt and solidification events to an ExaCA cellular-automata calculation. Define the sampling region with `box`, the Cartesian spacing with `dx`, and the crossing temperature with `isoValue`. The object interpolates AdditiveFOAM temperatures onto that lattice, so the number of sample points—and therefore memory and runtime—scales with the box volume divided by `dx³`. Reconstruct the per-process CSV files with `reconstructExaCAData` before supplying the result to ExaCA.

```foam
ExaCA
{
    libs    ("libadditiveFoamFunctionObjects.so");
    type    ExaCA;
    enabled true;
    box     (0 -1e-4 -2e-4) (0.002 1e-4 0);
    dx      2.5e-6;
    isoValue 1620;
}
```

`box` and sampling spacing `dx` define a Cartesian sampling lattice, and `isoValue` $$T_I$$ defaults to the liquidus. Temperatures at a sampling point are interpolated from its containing finite-volume cell. For a stored interval $$[t_0,t_1]$$ with endpoint temperatures $$T_0,T_1$$ that bracket $$T_I$$, the crossing fraction and time are

$$m=\operatorname{clip}_{[0,1]}\!\left[\frac{T_I-T_0}{T_1-T_0}\right],
\qquad t_I=t_0+m(t_1-t_0).$$

Here $$T_0$$ and $$T_1$$ are the interpolated point temperatures at times $$t_0$$ and $$t_1$$, $$m$$ is the fractional crossing position within that interval, and $$t_I$$ is the interpolated isovalue-crossing time.

An upward crossing stores the melt time $$t_m=t_I$$. The next downward crossing stores the solidification time $$t_s=t_I$$ and cooling rate

$$R=\frac{T_0-T_1}{t_1-t_0}.$$

Each rank writes `ExaCA/data_<rank>.csv` with columns `x,y,z,tm,ts,cr`. `reconstructExaCAData` creates `ExaCA/time-temperature.csv`.

{: .warning }
ExaCA supports mesh motion and topology-changing AMR on a fixed decomposition, but it deliberately aborts on dynamic mesh redistribution. Disable the ExaCA function object when using runtime load balancing.

The ExaCA sampling cost and memory use grow with the box volume divided by `dx³`.
