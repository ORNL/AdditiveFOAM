---
layout: versioned
title: Adaptive Mesh Refinement
parent: User Guide
nav_order: 7
permalink: /docs/amr/
redirect_from:
  - /reference/amr/
usemathjax: true
---

# Adaptive Mesh Refinement

AdditiveFOAM generates `refinementField` from two criteria: temperature above `refinementTemperature` and intersection with buffered, positive-power portions of the future scan path. OpenFOAM's `refiner` uses this field to refine and coarsen the mesh; its `loadBalancer` can then redistribute the changed mesh between MPI processes.

The selected refinement model determines the future scan-path interval represented by `refinementField`:

| Refinement model | Scan-path interval |
|---|---|
| **`timeStep`** | The next CFD time step |
| **`uniformTimeIntervals`** | The next interval in a fixed partition of the scan duration |
| **`targetCellLoad`** | A variable look-ahead interval chosen to approach a target number of cells per process |
| **`none`** | No AdditiveFOAM refinement marker |

The `refinementModel` dictionary is optional. Users select the model and source-specific buffer dimensions in `constant/heatSourceDict`, then configure refinement limits and optional load balancing in `constant/dynamicMeshDict`.

To enable AMR:

1. Select `timeStep`, `uniformTimeIntervals`, or `targetCellLoad` in `constant/heatSourceDict` and provide one buffer for each named heat source.
2. Configure OpenFOAM's `refiner` in `constant/dynamicMeshDict`, including `refineInterval 1`, `maxRefinement`, and `maxCells`.
3. Add the `loadBalancer` and a redistribution method in `system/decomposeParDict` when the refined mesh should be redistributed in parallel.
4. During the run, inspect the reported cell count and load imbalance. With `dumpLevel true`, reconstruct and visualize `cellLevel` and `refinementField` in ParaView to confirm where refinement is requested and applied.

<figure class="documentation-figure">
  <video controls autoplay loop muted playsinline poster="{{ '/assets/images/visualizations/amr-refinement.png' | relative_url }}" aria-label="Animation of adaptive mesh refinement following the moving heat source in the AMB2018-02-B case.">
    <source src="{{ '/assets/images/visualizations/amr-refinement.mp4' | relative_url }}" type="video/mp4">
    <img src="{{ '/assets/images/visualizations/amr-refinement.png' | relative_url }}" alt="Top view of the AMB2018-02-B mesh with the requested refinement region shown in green and the refined cells visible around it.">
  </video>
  <figcaption>Dynamic AMR in the AMB2018-02-B case. Green identifies cells selected by <code>refinementField</code>; the mesh color and edges show the resulting <code>cellLevel</code>. The OpenFOAM refiner is evaluated every time step with <code>refineInterval 1</code>.</figcaption>
</figure>

Active models require:

```foam
refinementModel
{
    refinementModel       timeStep;
    refinementTemperature 1000;
    buffers
    {
        beam (85e-6 85e-6 100e-6);
    }
}
```

| Entry | Required | Default | Meaning |
|---|---:|---|---|
| `refinementModel` | Yes in dictionary | — | `none`, `timeStep`, `uniformTimeIntervals`, or `targetCellLoad` |
| `refinementTemperature` | No | Disabled (`GREAT`) | Set $$R_c^T=1$$ where $$T_c\geq T_{\mathrm{ref}}$$ |
| `buffers` | Yes for active models | — | One scan-path buffer per source |

`maxRefinement` is read from `constant/dynamicMeshDict` and is required by active models.

## Refinement criterion

Let $$R_c\in\{0,1\}$$ be `refinementField` in cell $$c$$. At each refinement update,

$$R_c=\max\left(R_c^T,R_c^S\right),$$

where the temperature indicator is

$$R_c^T=H(T_c-T_{\mathrm{ref}})$$

and $$H$$ is the Heaviside function with $$H(0)=1$$. The scan-path indicator $$R_c^S$$ is one when the cell overlaps at least one buffered, positive-power path segment in the selected future time interval.

For a path segment from $$\mathbf{p}_0$$ to $$\mathbf{p}_1$$, define the horizontal scan frame

$$\mathbf{e}_0=\frac{(p_{1x}-p_{0x},p_{1y}-p_{0y},0)}
{\sqrt{(p_{1x}-p_{0x})^2+(p_{1y}-p_{0y})^2}},
\quad
\mathbf{e}_1=(-e_{0y},e_{0x},0),
\quad
\mathbf{e}_2=(0,0,1).$$

The oriented buffer has centre $$\mathbf{m}=(\mathbf{p}_0+\mathbf{p}_1)/2$$ and half-lengths

$$\mathbf{L}=
\left(\frac{\|\mathbf{p}_1-\mathbf{p}_0\|_{xy}}{2}+b_0,
b_1,b_2\right),$$

where `(b0 b1 b2)` is `(endpointPadding transverseHalfSpan verticalHalfSpan)`. For cell axis-aligned bounding-box half-lengths $$\mathbf{h}_c$$, its projected radius on buffer axis $$k$$ is

$$r_{c,k}=\sum_{j=1}^{3}h_{c,j}|e_{k,j}|.$$

Here $$\mathbf{p}_0$$ and $$\mathbf{p}_1$$ are the segment endpoints, $$\mathbf{e}_0$$ is the horizontal segment direction, $$\mathbf{e}_1$$ is its in-plane normal, and $$\mathbf{e}_2$$ is the vertical direction. The vector $$\mathbf{m}$$ is the segment midpoint, $$L_k$$ is the buffer half-length along axis $$k$$, $$h_{c,j}$$ is the cell bounding-box half-length along Cartesian direction $$j$$, and $$r_{c,k}$$ is that box's projected half-width on axis $$k$$. The notation $$\|\cdot\|_{xy}$$ denotes horizontal Euclidean length.

The overlap condition used to set $$R_c^S=1$$ is

$$\left|(\mathbf{c}_c-\mathbf{m})\cdot\mathbf{e}_k\right|
\le L_k+r_{c,k},\qquad k=0,1,2,$$

where $$\mathbf{c}_c$$ is the centre of the cell bounding box. Only path intervals with $$P_b>0$$ enter this union. The OpenFOAM `refiner` consumes $$R_c$$; mesh redistribution is a separate operation.

In the indicator equations, $$T_c$$ is the cell temperature, $$T_{\mathrm{ref}}$$ is `refinementTemperature`, $$R_c^T$$ is the temperature indicator, $$R_c^S$$ is the scan-path indicator, and $$P_b$$ is the power of source $$b$$.

## `targetCellLoad`

Select `targetCellLoad` when the future scan-path buffer should expand or contract in response to a desired mesh size per MPI process. Set `targetCellsPerProc` to the intended average workload and use `nBufferVolumes` to retain a minimum amount of upcoming path. The model adjusts its look-ahead time after each refinement update; `maxTargetVolumeGrowth` and `maxTargetVolumeShrink` limit how quickly that marked volume changes.

In `constant/heatSourceDict`:

```foam
refinementModel
{
    refinementModel       targetCellLoad;
    refinementTemperature 1000;

    buffers
    {
        beam (85e-6 85e-6 100e-6);
    }

    targetCellLoadCoeffs
    {
        targetCellsPerProc        5000;
        nBufferVolumes            4;
        maxSearchIter             30;
        timeTolerance             1e-6;
        initialTargetVolumeFactor 0.5;
        maxTargetVolumeGrowth     1.2;
        maxTargetVolumeShrink     0.8;
        postScanUpdateInterval    10;
    }
}
```

Let $$N$$ be the global cell count, $$N_p$$ the number of MPI processes, $$\overline V=V_\Omega/N$$ the mean cell volume, and $$L_{\max}$$ the maximum refinement level. Refining one cell through all levels adds

$$M=\max\left(2^{3L_{\max}}-1,1\right)$$

cells. For source-buffer half-lengths $$\mathbf{b}_s$$, the minimum projected volume is

$$V_{\min}=n_{\mathrm{buf}}\sum_s8b_{s,0}b_{s,1}b_{s,2}.$$

Here $$V_\Omega$$ is the mesh-domain volume, $$M$$ is the maximum cell-count increase per selected cell, $$s$$ indexes the heat sources, $$(b_{s,0},b_{s,1},b_{s,2})$$ are their buffer half-lengths, and $$n_{\mathrm{buf}}$$ is `nBufferVolumes`.

The minimum corresponding cell increment is $$N_{\min}=MV_{\min}/\overline V$$. Consequently, the requested cells per process are replaced by

$$N_{\mathrm{pp}}^{\star}=\max\left(N_{\mathrm{pp}},
\frac{N+N_{\min}}{N_p}\right).$$

Here $$N_{\mathrm{pp}}$$ is `targetCellsPerProc`, $$N_{\mathrm{pp}}^{\star}$$ is its minimum-volume-adjusted value, and $$N^{\star}=N_pN_{\mathrm{pp}}^{\star}$$ is the corresponding global target. With initial factor $$f_0$$ equal to `initialTargetVolumeFactor`, the initial target refinement volume is

$$V_R^{(0)}=\max\left[
f_0\overline V\frac{\max(N^{\star}-N,0)}{M},V_{\min}\right].$$

At update $$m$$, define the cell-count ratio

$$r^{(m)}=\frac{N_{\mathrm{pp}}^{\star}}{N^{(m)}/N_p}.$$

The target volume update is

$$V_R^{(m+1)}=\max\left[
V_R^{(m)}\operatorname{clip}_{[f_{\mathrm{shrink}},\,f_{\mathrm{grow}}]}
\!\left[r^{(m)}\right],
V_{\min}\right].$$

Here $$N^{(m)}$$ is the observed global cell count, $$r^{(m)}$$ is the target-to-observed cell-count ratio, and $$V_R^{(m)}$$ is the target marked volume at refinement update $$m$$. The quantities $$f_{\mathrm{shrink}}$$ and $$f_{\mathrm{grow}}$$ are `maxTargetVolumeShrink` and `maxTargetVolumeGrowth`, respectively.

Starting at the current time $$t$$, the look-ahead endpoint $$t_R$$ is the earliest positive-power path-event time for which

$$\sum_cV_cR_c(t,t_R)\ge V_R^{(m+1)}$$

and the newly marked scan-path volume is at least $$V_{\min}$$. Event intervals are accumulated in time order; the final interval is resolved by bisection to `timeTolerance` or `maxSearchIter`.

Here $$V_c$$ is cell volume, $$R_c(t,t_R)$$ is the marker evaluated over the future interval $$[t,t_R]$$, $$t$$ is the current time, and $$t_R$$ is the selected look-ahead endpoint.

## Mesh refinement and load balancing

In `constant/dynamicMeshDict`:

```foam
topoChanger
{
    type              refiner;
    libs              ("libfvMeshTopoChangers.so");
    refineInterval    1;
    field             refinementField;
    lowerRefineLevel  0.1;
    upperRefineLevel  10;
    nBufferLayers     5;
    maxRefinement     1;
    maxCells          20000000;
    dumpLevel         true;
}

distributor
{
    type                    loadBalancer;
    libs                    ("libfvMeshDistributors.so");
    redistributionInterval  10;
    maxImbalance            0.10;
}
```

`refineInterval` must be `1` for the AdditiveFOAM refinement models. The OpenFOAM `refiner` therefore evaluates the updated `refinementField` every CFD time step. `maxRefinement` supplies $$L_{\max}$$ in the equations above. `maxCells` bounds the global mesh size, and `redistributionInterval` is measured in time steps.

The load balancer reads its distributor from `system/decomposeParDict`:

```foam
numberOfSubdomains 8;

method scotch;

distributor zoltan;
libs ("libzoltanRenumber.so");

zoltanCoeffs
{
    lb_method   graph;
    lb_approach repartition;
}
```

For processor CPU-load estimates $$C_p$$, OpenFOAM defines the average $$\overline C=N_p^{-1}\sum_pC_p$$ and a conservative maximum $$C_{\max}$$ from the maximum base-CFD and registered load components. Redistribution occurs when

$$I=\frac{C_{\max}-\overline C}{\overline C}
>I_{\max},$$

where $$I_{\max}$$ is `maxImbalance`. Zoltan operates only in parallel and requires `libzoltanRenumber.so`; `lb_approach repartition` computes a new weighted decomposition from the current one.

Here $$p$$ indexes MPI processes, $$C_p$$ is the estimated load on process $$p$$, $$\overline C$$ is the process-average load, $$C_{\max}$$ is the conservative maximum load, and $$I$$ is the relative load imbalance.

## `timeStep`

Select `timeStep` when refinement should cover only the positive-power path traversed during the next CFD step. This gives the shortest look-ahead interval and has no model-specific coefficient dictionary.

At time $$t^n$$, this model evaluates the marker over $$[t^n,t^n+\Delta t]$$:

$$R_c=\max\!\left[H(T_c-T_{\mathrm{ref}}),R_c^S(t^n,t^n+\Delta t)\right].$$

Here $$R_c^S(t_0,t_1)=1$$ if cell $$c$$ intersects the oriented buffer of any positive-power scan-path segment active in $$[t_0,t_1]$$, and is zero otherwise. It has no model-specific coefficient dictionary. The segment–cell intersection test is defined in [Refinement criterion](#refinement-criterion).

## `uniformTimeIntervals`

Select `uniformTimeIntervals` when the scan should be divided into a fixed number of equal refinement windows. A smaller `intervals` value produces a longer look-ahead window and a larger simultaneously refined path volume; a larger value produces shorter windows. Configure the required count with:

```foam
uniformTimeIntervalsCoeffs
{
    intervals 100;
}
```

For scan end time $$t_e$$, start time $$t_s$$, and `intervals` $$N_I$$, the marker interval is

$$\Delta t_R=\frac{t_e-t_s}{N_I}.$$

Here $$\Delta t_R$$ is the uniform refinement-update duration and $$N_I$$ is the configured number of intervals.

At the $$j$$th update it evaluates $$R_c^S(t_j,t_j+\Delta t_R)$$. After $$t_e$$, $$R_c^S=0$$ and the scheduled update contains only $$H(T_c-T_{\mathrm{ref}})$$, permitting cells below the temperature threshold to be selected for coarsening.

## `none`

`none` disables calculation of $$R_c$$.

{: .warning }
The `ExaCA` function object aborts if it receives a dynamic mesh redistribution event. It supports refinement on a fixed processor decomposition, but it cannot currently be combined with runtime load balancing.
