---
layout: versioned
title: Projected Heat-source Calibration
parent: User Guide
nav_order: 5
permalink: /docs/heat-source-calibration/
usemathjax: true
---

# Projected Heat-source Calibration

`calibrateHeatSource` estimates the two coefficients that relate melt-pool depth to the axial shape of AdditiveFOAM's projected heat sources. It runs a collection of AdditiveFOAM cases, builds a local AdditiveFOAM response curve for every experimental condition, infers one local shape parameter per condition, and then fits the global projected-source closure. The command also caches completed simulations and writes a PDF report, a machine-readable fit, and a CSV summary.

This workflow calibrates the common depth projection used by `projectedGaussian`, `nLightAFX`, and `tabulated`. It does not calibrate the lateral beam profile, source dimensions, absorptivity, material properties, or boundary conditions. Those inputs must be established before interpreting the fitted depth-distribution coefficients.

The [heat-source model reference]({{ '/docs/heat-sources/#common-projected-source-formulation' | relative_url }}) defines the complete normalized volumetric distribution. This section focuses on the inverse problem solved by the calibration command.

## Model being calibrated

The projected models use the axial function

$$p(z;k,d_z)=\exp\left[-3\left(\frac{z}{d_z}\right)^k\right],$$

with normalization

$$\mathcal A_p(k,d_z)=\int_0^\infty p(z;k,d_z)\,\mathrm dz
=\frac{d_z\Gamma(1/k)}{k3^{1/k}}.$$

Here $$z$$ is distance into the material, $$d_z$$ is the current source depth, and $$k$$ controls the depth-profile shape. AdditiveFOAM computes $$k$$ through

$$x_s=\max\left[\frac{d_z}{s_0},10^{-3}\right],$$

$$n=\operatorname{clip}_{[0,9]}\!\left[A\log_2(x_s)+B\right],
\qquad k=2^n,$$

where

$$s_0=\min(d_x^0,d_y^0)$$

is the smaller configured lateral source dimension. The superscript $$0$$ identifies the dimensions read when the heat-source model is constructed; transient depth updates change $$d_z$$ but not $$s_0$$.

Increasing $$n$$ increases $$k$$, which makes the source more nearly uniform through most of its depth and sharpens its decay near $$d_z$$. The clipping limits imply $$1\leq k\leq512$$.

## Calibration overview

The calculation has two stages:

1. For each experimental process condition, set $$A=0$$ and infer a local value $$n_c=B_c$$ by comparing measured melt-pool depths with AdditiveFOAM trials.
2. Regress the local estimates against normalized measured depth to obtain the production coefficients $$A$$ and $$B$$.

The index $$c$$ denotes an experimental condition, such as one combination of power, scan speed, and spot size. The index $$r$$ denotes repeated depth measurements at that condition, and $$j$$ denotes a trial value used in an AdditiveFOAM response curve.

## Experimental observations

Each entry in `experiments.yml` supplies a process-condition dictionary and one or more measured depths:

```yaml
- parameters:
    Power_W: 300.0
    Speed_mm_s: 500
    Spot_Size_microns: 109.69
  Measured_Depth_microns: [168.57, 190.67]
```

`Spot_Size_microns` is the measured full D4σ beam diameter. For condition $$c$$, denote it by $$D_c$$. The command validates that it is finite and positive, converts it to the 2σ radius

$$w_c=\frac{D_c}{2},$$

and renders $$w_c$$ in metres through the `<<heatSourceRadius>>` placeholder. The supplied template assigns it to both lateral entries in `dimensions`, ensuring that every trial uses the beam size associated with its experimental condition.

The statistical calibration normalizes measured depth $$d_{c,r}$$ by that same lateral heat-source dimension:

$$y_{c,r}=\frac{d_{c,r}}{w_c}
=\frac{d_{c,r}}{D_c/2}.$$

## Stage 1: local AdditiveFOAM calibration

### Trial simulations

At condition $$c$$, every AdditiveFOAM trial fixes $$A=0$$ and sets the configured calibration token to a trial value $$b_j$$. The projected-source relation therefore reduces to

$$n=b_j,\qquad k=2^{b_j}.$$

The command copies the case template, replaces the placeholders for trial value, 2σ source radius, power, velocity, end time, and write interval, and runs the configured case command. Spot size is converted from micrometres to metres, and scan speed is converted from mm/s to m/s.

The `meltPoolDimensions` Function Object writes the isotherm history. The response recorded for trial $$j$$ is the maximum penetration over the simulated track,

$$\widehat d_{c,j}=\max_t d_{\mathrm{iso},c,j}(t),$$

and its normalized value is

$$\widehat y_{c,j}=\frac{\widehat d_{c,j}}{w_c}
=\frac{\widehat d_{c,j}}{D_c/2}.$$

When `case.melt_pool_isovalue` is `liquidus` or `solidus`, `foamDictionary` expands `constant/transportProperties`, reads `thermoPath`, and selects the temperature paired with `alpha.solid = 0` or `alpha.solid = 1`, respectively. A numeric temperature may be supplied instead.

### Response surrogate

The sampled pairs $$(b_j,\widehat y_{c,j})$$ are sorted by $$b_j$$. With at least three unique trials, the command evaluates a shape-preserving piecewise-cubic Hermite interpolant on a dense grid:

$$\widehat y_c(n)=\operatorname{PCHIP}
\left(\{b_j,\widehat y_{c,j}\}\right).$$

With only two unique trials, the sampled points are used directly and the likelihood interpolates linearly between them. The surrogate does not extrapolate beyond the smallest and largest completed trial values. PCHIP preserves local curve shape, but it does not force a non-monotonic AdditiveFOAM response to become monotonic.

### Bayesian inverse problem

The local prior is uniform over the simulated range:

$$n_c\sim\mathcal U\!\left(b_{\min},b_{\max}\right).$$

Repeated normalized depth measurements are represented with independent normal likelihoods:

$$y_{c,r}\mid n_c\sim
\mathcal N\!\left(\widehat y_c(n_c),\sigma_c^2\right),$$

where

$$\sigma_c=\max\left[
0.05\,\overline y_c,
\operatorname{std}\!\left(\{y_{c,r}\}\right)
\right].$$

The five-percent floor prevents identical or nearly identical repeated measurements from producing a zero-width likelihood. PyMC samples the posterior using the configured draws, tuning steps, cores, target acceptance, and random seed.

The reported local estimate is the retained posterior sample with the largest sampled log probability. The state file also records the mean, median, variance, standard deviation, 2.5th and 97.5th percentiles, and up to 1,000 posterior samples.

<figure class="documentation-figure documentation-figure--plot">
  <img src="{{ '/assets/images/visualizations/heat-source-calibration-responses.png' | relative_url }}" alt="Five projected heat-source calibration response curves showing AdditiveFOAM trial depth against local shape parameter, measured replicate ranges, and inferred local posterior modes.">
  <figcaption>Local calibration for the five tutorial powers. Circles are completed AdditiveFOAM trials, solid curves are the PCHIP surrogates, purple bands span the repeated measurements, and stars identify the inferred local posterior modes. The low- and high-<em>n</em> plateaus make the finite trial range and possible non-monotonic response visible.</figcaption>
</figure>

### Optional response-curve refinement

After the first posterior is sampled, the driver may add an AdditiveFOAM point when all of the following remain true:

- Fewer than `max_simulations_per_experiment` trials exist.
- The posterior standard deviation exceeds `posterior_std_tolerance`.
- No completed trial matches the mean measured depth within `depth_tolerance_microns`.
- The target depth crosses the PCHIP response curve within `calibration.bounds`.

The new trial is obtained by interpolating the first response-curve crossing of the mean measured depth. The AdditiveFOAM case is run, the surrogate is rebuilt, and the posterior is sampled again. Refinement stops as soon as either uncertainty or depth error satisfies its tolerance, the trial limit is reached, or no in-range crossing exists.

If the initial-value list already contains `max_simulations_per_experiment` values, as in the supplied tutorial, the response curve is an exhaustive fixed grid and no additional trial can be added.

## Stage 2: global closure fit

For every completed condition, the global fit uses the mean normalized measured depth

$$x_c=\frac{1}{N_c}\sum_{r=1}^{N_c}\frac{d_{c,r}}{w_c},
\qquad z_c=\log_2(x_c),$$

the local posterior mode $$n_c$$, and its posterior standard deviation $$s_{n,c}$$. The standardized residual is

$$r_c(A,B)=\frac{n_c-(Az_c+B)}{s_{n,c}}.$$

Weighted linear least squares supplies the initial estimate. With at least three conditions, SciPy then minimizes the robust objective

$$\frac12\sum_c\rho\!\left(r_c^2\right),$$

using `soft_l1` by default. The uncertainty weighting gives better-constrained local calibrations more influence, while the robust loss limits the influence of an isolated condition that is inconsistent with the global trend.

The approximate parameter covariance is

$$\operatorname{Cov}(A,B)\approx
(J^\mathsf TJ)^{-1}\max(\chi_r^2,1),$$

where $$J$$ is the uncertainty-scaled design matrix and $$\chi_r^2$$ is the reduced weighted residual sum of squares.

The reported production relation is

$$n(x)=\operatorname{clip}_{[0,9]}\!\left[A\log_2(x)+B\right],
\qquad k(x)=2^{n(x)}.$$

The robust fit itself is performed on the unclipped linear relation. Clipping is applied when the response curve is evaluated and when AdditiveFOAM uses the coefficients.

### Empirical uncertainty band

The final response interval propagates the local posterior uncertainty through complete refits:

1. Draw one $$n_c$$ from each stored local posterior. If samples are unavailable, draw from a normal approximation using the stored local mode and standard deviation.
2. Clip the sampled local values to $$[0,9]$$.
3. Refit $$A$$ and $$B$$ with the same robust, uncertainty-weighted procedure.
4. Evaluate and clip the response curve over the plotted $$x$$ range.

The pointwise 2.5th, 50th, and 97.5th percentiles of the successful bootstrap curves form the empirical interval. The report also records corresponding percentile intervals for $$A$$ and $$B$$. If too few bootstrap fits succeed, the plot falls back to the covariance approximation.

<figure class="documentation-figure documentation-figure--plot">
  <img src="{{ '/assets/images/visualizations/heat-source-calibration-fit.png' | relative_url }}" alt="Final projected heat-source calibration relation showing five local estimates, asymmetric uncertainty bars, the robust logarithmic fit, and an empirical 95 percent interval.">
  <figcaption>The tutorial's final 2σ-normalized closure, with x = d/(D4σ/2). Points and error bars are local posterior modes and 95% intervals; the green curve is the robust uncertainty-weighted fit, and the orange band propagates local posterior uncertainty through repeated global refits.</figcaption>
</figure>

## Normalization contract

`Spot_Size_microns` is the measured full D4σ diameter. The driver computes the 2σ heat-source radius

$$w=\frac{D_{4\sigma}}{2},$$

renders it into both lateral source dimensions, and uses it in the calibration coordinate:

$$x=\frac{d_z}{w}
=\frac{d_z}{D_{4\sigma}/2}.$$

The projected heat-source models evaluate the same coordinate:

$$x=\frac{d_z}{\min(d_x^0,d_y^0)}.$$

For the supplied circular source, $$d_x^0=d_y^0=w$$, so the calibration and runtime relations are identical.

{: .important }
Use `A` and `B` from `calibration_fit.yml` directly as the source slope and source intercept.

## Configuration

`config.yml` contains five functional groups:

| Section | Entry | Meaning |
|---|---|---|
| `paths` | `experiments` | Experimental YAML file |
| `paths` | `template` | Complete AdditiveFOAM case containing renderer placeholders |
| `paths` | `campaign` | Directory for generated cases, caches, fit, and reports |
| `case` | `command` | Command run inside each rendered case |
| `case` | `keep_successful` | Retain or remove processor and numeric time directories after metrics are saved |
| `case` | `melt_pool_isovalue` | `liquidus`, `solidus`, or a temperature in kelvin |
| `calibration` | `token` | Placeholder name used for the local trial parameter |
| `calibration` | `initial_values` | Initial local values used to form each response curve |
| `calibration` | `bounds` | Permitted local-value interval |
| `calibration` | `max_simulations_per_experiment` | Maximum AdditiveFOAM trials for one condition |
| `calibration` | `depth_tolerance_microns` | Stop refinement when a trial is this close to mean measured depth |
| `calibration` | `posterior_std_tolerance` | Stop refinement when local posterior uncertainty falls below this value |
| `calibration` | `pchip_grid_points` | Dense-grid resolution for the response surrogate |
| `bayesian` | `draws`, `tune`, `cores` | PyMC sampling controls |
| `bayesian` | `target_accept` | Sampler target acceptance probability |
| `bayesian` | `random_seed` | Reproducible local-posterior seed |
| `final_fit` | `bootstrap_samples` | Number of uncertainty-propagation refits |
| `final_fit` | `bootstrap_random_seed` | Reproducible bootstrap seed |
| `report` | `enabled` | Generate the PDF and CSV report products |

Relative paths are resolved from the directory containing `config.yml`. Environment variables and `~` are expanded.

## Case-template contract

The template must contain `0/`, `constant/heatSourceDict`, `constant/scanPath`, `constant/transportProperties`, `system/controlDict`, and `Allrun`. The renderer requires these placeholders exactly once:

| File | Placeholder |
|---|---|
| `constant/heatSourceDict` | The token selected by `calibration.token`, such as `<<B>>`, and `<<heatSourceRadius>>` |
| `constant/scanPath` | `<<power>>` and `<<velocity>>` |
| `system/controlDict` | `<<endTime>>` and `<<writeInterval>>` |

The current end-time calculation sums the travel time of mode-`0` linear segments. It does not add mode-`1` dwell durations. Calibration templates should therefore use the supplied initial-position row followed by translating segments, unless `endTime` handling is extended for dwell-based paths.

## Cache and output files

The campaign directory is persistent:

```text
campaign/
├── cases/<condition>/<trial>/
│   ├── log.run
│   └── metrics.yml
├── simulations.yml
├── calibration_state.yml
├── calibration_fit.yml
└── reports/
    ├── calibration_report.pdf
    └── calibration_summary.csv
```

`simulations.yml` stores completed trial values and simulated depths by process-condition dictionary. `calibration_state.yml` stores local posteriors; a condition is considered current when its parameter dictionary matches and the fingerprint of its measured-depth list is unchanged. `calibration_fit.yml` stores the final coefficients, diagnostics, covariance approximation, bootstrap intervals, and calibrated $$x$$ range.

{: .warning }
The cache fingerprint does not include the case template, material files, solver revision, mesh, calibration settings, or Bayesian settings. After changing any of those inputs, remove the affected cached simulation and state data, or choose a new campaign directory, before rerunning the calibration.

With `keep_successful: false`, rendered inputs, `log.run`, Function Object output, and `metrics.yml` remain, but processor and numeric time directories are removed after the maximum melt-pool depth is recorded.

The complete worked example is the [heat-source calibration tutorial]({{ '/tutorials/heat-source-calibration/' | relative_url }}).
