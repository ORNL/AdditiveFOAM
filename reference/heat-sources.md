---
layout: versioned
title: Heat Source Models
parent: User Guide
nav_order: 4
permalink: /docs/heat-sources/
redirect_from:
  - /reference/heat-sources/
usemathjax: true
---

# Heat Source Models

AdditiveFOAM represents deposited energy with moving volumetric heat sources. Each source follows an arbitrary piecewise-linear or dwell-based scan path that specifies position, time, and power. A source independently selects an absorption model, which determines the absorbed fraction of the prescribed power, and a spatial distribution, which determines how that absorbed power is distributed through the material. Multiple sources can operate in the same simulation with separate paths, powers, absorption models, dimensions, and integration controls. The resulting volumetric power density is normalized so its spatial integral equals the absorbed source power.

## Scan paths

The `pathName` entry selects a file under `constant/`. Its first line is a header and is ignored. Each subsequent non-empty row has six columns:

```text
mode x y z power parameter
```

| Mode | Position columns | Power column | Parameter column |
|---|---|---|---|
| `0` | Endpoint of a linear segment, in metres | Source power during the segment, in watts | Translation speed, in m/s |
| `1` | Fixed source position, in metres | Source power during the dwell, in watts | Dwell duration, in seconds |

For mode `0`, the source moves from the preceding path position to $$(x,y,z)$$ at the specified speed. For mode `1`, it remains at $$(x,y,z)$$ for the specified duration. A zero-power row represents motion or dwell without energy deposition. The path begins at simulation time zero; an initial mode-`1` row can establish the starting position and an optional delay before the first translating segment.

For example,

```text
Mode  X       Y        Z    Power  Parameter
1     0.000   0.0000   0.0  0.0    0.0
0     0.002   0.0000   0.0  195.0  0.8
0     0.002   0.0001   0.0  0.0    0.8
0     0.000   0.0001   0.0  195.0  0.8
```

defines two powered, 2 mm scan segments at 0.8 m/s separated by a zero-power 0.1 mm transverse move. The [scan-path file section]({{ '/docs/case-files/#scan-path-file' | relative_url }}) gives the corresponding time interpolation equations.

## Top-level structure

Create `constant/heatSourceDict` and list every source name in `sources`. Each named dictionary must select its scan-path file, absorption model, and spatial distribution. Use `deltaT` to control temporal sampling of source motion and `nPoints` inside the distribution coefficients to control spatial quadrature. The example below defines one source named `beam1`; additional sources use the same structure with unique names and paths.

```foam
sources (beam1 beam2);

beam1
{
    pathName         scanPath_1;
    deltaT           1e-7;
    hitPathIntervals true;

    absorptionModel  constant;
    constantCoeffs { eta 0.33; }

    heatSourceModel  superGaussian;
    superGaussianCoeffs
    {
        dimensions (85e-6 85e-6 30e-6);
        nPoints    (10 10 10);
        k          2;
    }
}
```

| Entry | Required | Default | Meaning |
|---|---:|---|---|
| `sources` | Yes | — | Ordered list of beam dictionary names |
| `pathName` | Yes | — | File under `constant/` |
| `deltaT` | No | Unlimited | Internal temporal sampling interval for beam motion |
| `hitPathIntervals` | No | `true` | Make the CFD step land on scan-path boundaries |
| `absorptionModel` | Yes | — | `constant` or `Kelly` |
| `heatSourceModel` | Yes | — | Spatial distribution type |
| `dimensions` | Yes | — | Characteristic `(x y z)` dimensions in metres |
| `nPoints` | No | `(1 1 1)` | Target sub-cell integration density |
| `transient` | No | `false` | Update source depth from an isotherm |
| `isoValue` | No | Material liquidus when transient | Temperature used for transient depth |

See [adaptive heat-source integration](#adaptive-heat-source-integration) for how `nPoints` and `deltaT` differ.

## Moving volumetric sources

For $$N_b$$ heat sources, the volumetric power density is

$$\dot q(\mathbf{x},t)
=\sum_{b=1}^{N_b}
\frac{\eta_b(a_b)P_b(t)}{V_{b,\mathrm{norm}}}
w_b\!\left(\mathbf{x}-\mathbf{x}_b(t)\right),$$

where $$P_b(t)$$ and $$\mathbf{x}_b(t)$$ are obtained by piecewise interpolation of the scan path, $$\eta_b$$ is absorptivity, $$a_b=d_{b,z}/\min(d_{b,x},d_{b,y})$$ is the source aspect ratio, $$w_b$$ is the selected spatial kernel, and $$V_{b,\mathrm{norm}}$$ normalizes its integral.

Here $$\mathbf{x}$$ is spatial position, $$t$$ is time, $$N_b$$ is the number of sources, $$b$$ is the source index, $$P_b$$ is power, $$\mathbf{x}_b$$ is the source centre, and $$(d_{b,x},d_{b,y},d_{b,z})$$ are the source dimensions.

### Adaptive heat-source integration

The solver integrates each spatial distribution over finite-volume cells and averages its motion over every CFD time step. Users control these two approximations independently: `nPoints` sets the target spatial sampling density and the source-level `deltaT` sets the maximum temporal sub-interval. Increase `nPoints` until `plotPower` shows that the volume-integrated source has converged to the absorbed power. Reduce `deltaT` when the source moves a significant distance within one CFD step and confirm that the resulting temperature history no longer changes at the required resolution.

The finite-volume source for cell $$c$$ is based on the cell-average kernel weight

$$\overline{w}_{b,c}=\frac{1}{V_c}\int_{\Omega_c}
w_b\!\left(\mathbf{x}-\mathbf{x}_b\right)\,\mathrm{d}V.$$

Here $$\Omega_c$$ is the volume occupied by cell $$c$$, $$V_c=\int_{\Omega_c}\mathrm{d}V$$ is its volume, and $$\overline w_{b,c}$$ is the cell-average kernel value.

#### Spatial integration

Every heat-source coefficient dictionary accepts `nPoints`:

```foam
superGaussianCoeffs
{
    k          2.0;
    dimensions (85e-6 85e-6 30e-6);
    nPoints    (10 10 10);
}
```

For source dimension $$d_{b,j}$$ and configured count $$n_{b,j}$$, the target sampling distance in direction $$j$$ is

$$h_{b,j}=\frac{d_{b,j}}{n_{b,j}}.$$

For a hexahedral cell with axis-aligned bounding-box span $$\ell_{c,j}$$, the number and width of integration intervals are

$$N_{c,j}=\max\left(1,\left\lfloor\frac{\ell_{c,j}+\varepsilon}{h_{b,j}}\right\rfloor\right),
\qquad
\Delta x_{c,j}=\frac{\ell_{c,j}}{N_{c,j}}.$$

Here $$j\in\{x,y,z\}$$ denotes a Cartesian direction, $$h_{b,j}$$ is the target sample spacing, $$\ell_{c,j}$$ is the cell bounding-box span, $$N_{c,j}$$ is the resulting interval count, $$\Delta x_{c,j}$$ is the interval width, and $$\varepsilon$$ is the small positive tolerance used before integer truncation.

Midpoint quadrature gives

$$\overline{w}_{b,c}
\approx\frac{\Delta V_c}{V_c}
\sum_{i=1}^{N_{c,x}}
\sum_{j=1}^{N_{c,y}}
\sum_{k=1}^{N_{c,z}}
w_b(\mathbf{x}_{ijk}-\mathbf{x}_b),
\qquad
\Delta V_c=\prod_{j=1}^{3}\Delta x_{c,j}.$$

The indices $$i$$, $$j$$, and $$k$$ enumerate midpoint samples in the three Cartesian directions; $$\mathbf{x}_{ijk}$$ is a sample position and $$\Delta V_c$$ is its quadrature volume.

For non-hexahedral cells, $$\overline{w}_{b,c}=w_b(\mathbf{x}_c-\mathbf{x}_b)$$ is evaluated at the cell centre. The discrete kernel integral is

$$V_{b,h}=\sum_c V_c\overline{w}_{b,c}.$$

Let $$V_{b,0}$$ denote the normalization supplied by the heat-source model. The normalization used in the source term is

$$V_{b,\mathrm{norm}}=
\begin{cases}
V_{b,h}, & \left|1-V_{b,h}/V_{b,0}\right|<0.05,\\
V_{b,0}, & \text{otherwise}.
\end{cases}$$

The cell source is therefore

$$\dot q_{b,c}=\eta_bP_b\frac{\overline{w}_{b,c}}{V_{b,\mathrm{norm}}},$$

and its integrated power is $$\sum_cV_c\dot q_{b,c}$$. `nPoints` controls the quadrature resolution. Increase it until this integral converges to $$\eta_bP_b$$; `plotPower` plots the integrated value.

Sub-cell quadrature is implemented for `superGaussian`, `modifiedSuperGaussian`, `projectedGaussian`, `nLightAFX`, and `tabulated`.

#### Time integration

Let the CFD step be $$[t^n,t^{n+1}]$$ and the configured beam integration interval be $$\delta t_b$$. It is partitioned into $$K$$ sub-intervals with

$$\delta t_k=\min\left(\delta t_b,
t^{n+1}-t^n-\sum_{j=1}^{k-1}\delta t_j\right).$$

The time-averaged source is

$$\overline{\dot q}_{b,c}^{\,n+1}
=\frac{1}{\Delta t}
\sum_{k=1}^{K}\delta t_k\,
\dot q_{b,c}\!\left(t^n+\sum_{j=1}^{k}\delta t_j\right),
\qquad
\Delta t=\sum_{k=1}^{K}\delta t_k.$$

Here $$n$$ is the CFD time-level index, $$K$$ is the number of source sub-intervals, $$k$$ is the sub-interval index, $$\delta t_k$$ is its duration, and $$\Delta t=t^{n+1}-t^n$$ is the CFD-step duration.

The beam-level `deltaT` entry specifies $$\delta t_b$$:

```foam
beam
{
    pathName         scanPath;
    deltaT           2e-7;
    hitPathIntervals true;
    // absorption and heat-source entries ...
}
```

With `hitPathIntervals true`, the remaining time $$\tau$$ to the next scan-path event is divided into an integer number of CFD steps, so a step endpoint coincides with that event.

## Absorption models

The absorption model converts prescribed scan-path power into absorbed power. Select `constant` for a fixed absorptivity or `Kelly` when absorptivity should vary with the source depth-to-width aspect ratio. The selected value multiplies the normalized spatial distribution and is reported through the integrated source power in the solver log.

### `constant`

```foam
absorptionModel constant;
constantCoeffs { eta 0.33; }
```

`eta` is required.

### `Kelly`

The Kelly model represents the increase in effective absorptivity caused by repeated reflections within a conical or cylindrical depression. It uses a fixed absorptivity below a prescribed depth-to-width ratio and the geometry-dependent multiple-reflection relation above it.

```foam
absorptionModel Kelly;
KellyCoeffs
{
    geometry          cone;
    eta0              0.28;
    etaMin            0.35;
    aspectRatioSwitch 1.0;
}
```

Let $$a=d_z/\min(d_x,d_y)$$ be the source aspect ratio and $$a_s$$ be `aspectRatioSwitch`. The implemented absorptivity is

$$
\eta(a)=
\begin{cases}
\eta_{\min}, & a\leq a_s,\\[4pt]
\displaystyle
\eta_0\frac{1+(1-\eta_0)[G(a)-F(a)]}
{1-(1-\eta_0)[1-G(a)]}, & a>a_s.
\end{cases}
$$

For $$a>a_s$$, define $$\theta=\tan^{-1}(1/a)$$. The geometric functions are

$$
[F(a),G(a)]=
\begin{cases}
\left[
\dfrac{3\sin\theta-\sin(3\theta)}{4},
\dfrac{1}{1+\sqrt{1+a^2}}
\right], & \texttt{geometry cone},\\[12pt]
\left[
\dfrac{1-\cos(2\theta)}{2},
\dfrac{1}{2(1+a)}
\right], & \texttt{geometry cylinder}.
\end{cases}
$$

`eta0`, `etaMin`, and `geometry` are required; $$a_s$$ defaults to one.

In these equations, $$\eta_0$$ is the base optical absorptivity, $$\eta_{\min}$$ is the value used below the switch, $$\theta$$ is the depression half-angle, and $$F$$ and $$G$$ are the selected cone or cylinder geometry factors.

<figure class="documentation-figure documentation-figure--plot">
  <img src="{{ '/assets/images/visualizations/kelly-absorption.png' | relative_url }}" alt="Kelly effective absorptivity as a function of source aspect ratio for cone and cylinder geometries.">
  <figcaption>Kelly absorptivity for the parameters in the example dictionary. The model uses the fixed value $\eta_{\min}$ through the aspect-ratio switch $a_s$, after which the cone or cylinder multiple-reflection relation determines $\eta(a)$.</figcaption>
</figure>

## Heat-source models

The heat-source model determines the three-dimensional distribution of absorbed power. Select the model according to the available beam description:

| Model | Spatial input |
|---|---|
| `superGaussian` | One analytic three-dimensional ellipsoidal distribution |
| `modifiedSuperGaussian` | Analytic lateral super-Gaussian whose width contracts with depth |
| `projectedGaussian` | Analytic elliptical Gaussian intensity projected through depth |
| `nLightAFX` | Two characterized Gaussian-ring components projected through depth |
| `tabulated` | Measured or computed lateral intensity table projected through depth |

<figure class="documentation-figure">
  <img src="{{ '/assets/images/visualizations/heat-source-models.png' | relative_url }}" alt="Top-surface and centre-plane distributions for the five AdditiveFOAM heat-source models.">
  <figcaption>Heat-source distributions for the five AdditiveFOAM models. The <code>nLightAFX</code> and <code>tabulated</code> distributions use the Index 6 profiles supplied with the tutorials.</figcaption>
</figure>

All models use `dimensions` and accept `nPoints`. For models that accept `transient true`, AdditiveFOAM updates only the depth dimension from the position of `isoValue`; if `isoValue` is omitted, the material liquidus is used. This allows the volumetric distribution and aspect-ratio-dependent absorption or projection to follow the evolving penetration depth.

### `superGaussian`

The `superGaussian` model defines one ellipsoidal distribution whose shape is uniform in the scaled coordinate $$\sum_j(r_j/s_j)^2$$. The exponent $$k$$ controls the transition from a Gaussian profile at $$k=2$$ to a flatter central region with steeper decay for $$k>2$$.

For displacement $$\mathbf{r}=\mathbf{x}-\mathbf{x}_b$$ and

$$s_j=\frac{d_j}{2^{1/k}},$$

the kernel and analytic normalization are

$$w(\mathbf{r})=\exp\left\{-\left[
\sum_{j=1}^{3}\left(\frac{r_j}{s_j}\right)^2
\right]^{k/2}\right\},$$

$$V_0=\frac{2}{3}\pi s_xs_ys_z
\Gamma\!\left(1+\frac{3}{k}\right).$$

Here $$r_j$$ is component $$j$$ of the source-relative displacement, $$d_j$$ is the corresponding configured dimension, $$s_j$$ is its scaled length, $$V_0$$ is the analytic kernel integral over the active half-space, and $$\Gamma$$ is the Euler gamma function.

Thus `k=2` is Gaussian. `k` is required.

```foam
superGaussianCoeffs
{
    dimensions (85e-6 85e-6 30e-6);
    k          2;
    nPoints    (10 10 10);
}
```

### `modifiedSuperGaussian`

The `modifiedSuperGaussian` model applies a super-Gaussian profile in each plane normal to the source-depth direction. Its lateral scale decreases continuously with $$\lvert r_z\rvert$$ and becomes zero at $$\lvert r_z\rvert=d_z$$, producing a bounded three-dimensional distribution.

Let $$s_x=d_x/2^{1/k}$$, $$s_y=d_y/2^{1/k}$$, $$s_z=d_z$$, and $$z=\lvert r_z\rvert$$. For $$z<s_z$$, define

$$g(z)=\left[1-\left(\frac{z}{s_z}\right)^m\right]^{1/m}$$

and

$$w(\mathbf{r})=\exp\left\{-\left[
\left(\frac{r_x}{s_xg}\right)^2+
\left(\frac{r_y}{s_yg}\right)^2
\right]^{k/2}\right\};$$

$$w=0$$ for $$z\geq s_z$$. The analytic normalization is

$$V_0=\pi s_xs_ys_z\Gamma\!\left(1+\frac{2}{k}\right)
\frac{\Gamma(1+1/m)\Gamma(1+2/m)}{\Gamma(1+3/m)}.$$

Here $$g(z)$$ is the depth-dependent lateral scale multiplier, $$V_0$$ is the kernel integral, `k` controls the lateral exponent, and `m` controls the contraction with depth; both coefficients are required.

```foam
modifiedSuperGaussianCoeffs
{
    dimensions (40e-6 40e-6 30e-6);
    k          7.95;
    m          2.72;
    transient  true;
    nPoints    (10 10 10);
}
```

### Common projected-source formulation

The `projectedGaussian`, `tabulated`, and `nLightAFX` models construct a volumetric distribution from a lateral intensity and a depth projection. Let $$(x,y,z)$$ be coordinates relative to the source centre, with $$z\geq0$$ directed into the active material. For lateral intensity $$I(x,y)$$ and projection $$p(z)$$, define the normalized distribution

$$Q(x,y,z)=
\frac{I(x,y)p(z)}{\mathcal{A}_I\mathcal{A}_p},$$

where

$$\mathcal{A}_I=\int_{-\infty}^{\infty}\int_{-\infty}^{\infty}
I(x,y)\,\mathrm{d}x\,\mathrm{d}y,
\qquad
\mathcal{A}_p=\int_0^\infty p(z)\,\mathrm{d}z.$$

Thus $$Q$$ has units of inverse volume and satisfies

$$\int_0^\infty\int_{-\infty}^{\infty}\int_{-\infty}^{\infty}
Q(x,y,z)\,\mathrm{d}x\,\mathrm{d}y\,\mathrm{d}z=1.$$

The volumetric power density contributed by one source is

$$\dot q(x,y,z,t)=\eta P(t)Q(x,y,z),$$

where $$P(t)$$ is the scan-path power and $$\eta$$ is the absorptivity. All three models use the projection

$$p(z;k,d_z)=\exp\left[-3\left(\frac{z}{d_z}\right)^k\right],
\qquad
\mathcal{A}_p(k,d_z)=\frac{d_z\Gamma(1/k)}{k\,3^{1/k}},$$

where $$d_z$$ is the source depth and $$k$$ is its projection exponent. Define the clipped logarithmic parameter

$$a=\max\left[\frac{d_z}{\min(d_x^0,d_y^0)},10^{-3}\right],
\qquad
n=\operatorname{clip}_{[0,9]}\!\left[A\log_2(a)+B\right].$$

Here $$a$$ is the depth-to-width aspect ratio, $$d_x^0$$ and $$d_y^0$$ are the configured lateral dimensions, $$A$$ and $$B$$ are the model coefficients, and $$n$$ is the clipped base-two logarithm of the projection exponent. The exponent used in $$p(z;k,d_z)$$ is $$k=2^n$$. A transient depth changes $$d_z$$ and therefore updates $$a$$, $$n$$, $$k$$, and $$\mathcal{A}_p$$.

<figure class="documentation-figure documentation-figure--plot">
  <img src="{{ '/assets/images/visualizations/heat-source-projection.png' | relative_url }}" alt="Normalized common depth-projection function for exponents 1, 2, 4, and 8.">
  <figcaption>The common projected-source depth function for representative exponents. Increasing the exponent produces a more uniform distribution through most of the configured depth followed by a sharper decay near the specified source depth.</figcaption>
</figure>

### `projectedGaussian`

The `projectedGaussian` model uses an elliptical Gaussian lateral intensity with the common depth projection above:

$$I_G(x,y)=\exp\left[-2\left(
\frac{x^2}{d_x^2}+\frac{y^2}{d_y^2}\right)\right],
\qquad
\mathcal{A}_{I,G}=\frac{\pi d_xd_y}{2},$$

where $$d_x$$ and $$d_y$$ are the lateral source dimensions. Its normalized distribution is $$Q_G=I_Gp/(\mathcal{A}_{I,G}\mathcal{A}_p)$$.

```foam
projectedGaussianCoeffs
{
    dimensions (85e-6 85e-6 30e-6);
    A          0;
    B          1;
    transient  true;
    nPoints    (10 10 10);
}
```

### `nLightAFX`

The `nLightAFX` model represents a measured beam as the weighted sum of two normalized concentric Gaussian-ring components. Each component uses the common projection-function form, but its coefficients can produce a different depth exponent.

For component $$i\in\{0,1\}$$, let $$r_i$$ be its ring radius, $$\sigma_i$$ its radial standard deviation, and $$A_i$$ and $$B_i$$ its projection coefficients. Its clipped logarithmic parameter is

$$n_i=\operatorname{clip}_{[0,9]}\!\left[A_i\log_2(a)+B_i\right].$$

The corresponding projection exponent is $$k_i=2^{n_i}$$.

With radial coordinate $$r=\sqrt{x^2+y^2}$$, define the lateral intensity and depth projection

$$I_i(r)=\exp\left[-\frac12\left(\frac{r-r_i}{\sigma_i}\right)^2\right]
+\exp\left[-\frac12\left(\frac{r+r_i}{\sigma_i}\right)^2\right],$$

$$p_i(z)=p(z;k_i,d_z).$$

The volume integral of component $$i$$ is

$$\mathcal{A}_i=
\left(\int_{\mathbb R^2}I_i\,\mathrm{d}x\,\mathrm{d}y\right)
\left(\int_0^\infty p_i\,\mathrm{d}z\right)
=\frac{2\pi\sigma_i d_z\Gamma(1/k_i)}{k_i3^{1/k_i}}
\left[2\sigma_i e^{-r_i^2/(2\sigma_i^2)}
+\sqrt{2\pi}r_i\operatorname{erf}\left(\frac{r_i}{\sqrt2\sigma_i}\right)\right].$$

For outer-component fraction `alpha` $$=\alpha$$, the normalized distribution is

$$Q_{\mathrm{AFX}}(r,z)
=(1-\alpha)\frac{I_0(r)p_0(z)}{\mathcal{A}_0}
+\alpha\frac{I_1(r)p_1(z)}{\mathcal{A}_1},$$

where $$0\leq\alpha\leq1$$ and indices $$0$$ and $$1$$ denote the inner and outer components, respectively.

Each component requires `radius`, `sigma`, `A`, and `B`; `dimensions` is shared.

```foam
nLightAFXCoeffs
{
    dimensions (142.905e-6 142.905e-6 50e-6);
    alpha      0.902;
    inner { radius 33.50e-6; sigma 18.510e-6; A 0; B 1; }
    outer { radius 101.41e-6; sigma 16.3475e-6; A 0; B 1; }
    transient true;
    nPoints   (10 10 10);
}
```

The tutorial-supplied `nLightAFX.cfg` contains characterized modes `Index0` through `Index6` that can be expanded with OpenFOAM dictionary substitution.

### `tabulated`

The `tabulated` model supplies $$I(x,y)$$ as a measured or computed regular-grid table and combines it with the common depth projection.

Let $$f_{ij}=I(x_i,y_j)$$ be the intensity at grid node $$(x_i,y_j)$$, and let $$\Delta x$$ and $$\Delta y$$ be the uniform grid spacings. Inside cell $$(i,j)$$, define local coordinates $$\xi=(x-x_i)/\Delta x$$ and $$\upsilon=(y-y_j)/\Delta y$$. The lateral intensity is the bilinear interpolant

$$f(x,y)=(1-\xi)(1-\upsilon)f_{ij}
+\xi(1-\upsilon)f_{i+1,j}
+(1-\xi)\upsilon f_{i,j+1}
+\xi\upsilon f_{i+1,j+1}.$$

It is zero outside the nodal domain. For $$n_x$$ nodes in $$x$$ and $$n_y$$ nodes in $$y$$, its lateral normalization is

$$\mathcal{A}_{I,T}=\frac{\Delta x\Delta y}{4}
\sum_{j=0}^{n_y-2}\sum_{i=0}^{n_x-2}
\left(f_{ij}+f_{i+1,j}+f_{i,j+1}+f_{i+1,j+1}\right).$$

The normalized volumetric distribution is $$Q_T=f(x,y)p(z)/(\mathcal{A}_{I,T}\mathcal{A}_p)$$. Here $$i$$ and $$j$$ are the grid indices, while $$n_x$$ and $$n_y$$ are the node counts read from the first line of the table.

```foam
tabulatedCoeffs
{
    file       "beamProfile.txt";
    dimensions (250e-6 250e-6 50e-6);
    A          0;
    B          1;
    nPoints    (10 10 10);
}
```

`file` is relative to `constant/`. Its headerless format is:

```text
nx ny
x0 y0
dx dy
f00 f10 ...
f01 f11 ...
```

Coordinates and spacing are metres. Values are row-major with x varying fastest. Use `primesToAdditiveFoam` to normalize compatible PRIMES CSV exports.
