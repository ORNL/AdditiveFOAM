# Tabulated Heat Source Tutorial

This tutorial demonstrates the `tabulated` heat source model in AdditiveFOAM. The model reads a user-defined, uniformly spaced 2D beam intensity distribution from an ASCII file and applies the same projected axial distribution used by the `projectedGaussian` heat source.

The purpose of this tutorial is to show how arbitrary beam shapes can be tested without adding a new analytic heat source model for every profile.

## File structure

The important file for this tutorial are:

```text
constant/beamProfile.txt
```

Defines the 2D planar beam intensity profile used by the `tabulated` heat source.

## Heat source model

The tutorial uses:

```foam
heatSourceModel tabulated;
```

The corresponding coefficient dictionary is:

```foam
tabulatedCoeffs
{
    dimensions  (2.50e-4 2.50e-4 5.0e-5);
    transient   true;
    isoValue    1700;
    nPoints     (4 4 4);

    A           1.69;
    B          -0.12;

    file        "beamProfile.txt";
}
```

### Coefficients

`dimensions`

Sets the heat source dimensions used by the base moving heat source integration. For this tutorial, the lateral dimensions are set to cover the tabulated beam support, and the third component is the initial projected depth.

`A` and `B`

Define the projected axial shape closure:

```text
n = A*log2(x) + B
```

where `x` is the ratio between the current heat source depth and lateral heat source size. The implementation clamps this internal numerical exponent consistently with the `projectedGaussian` model.

`file`

Name of the tabulated 2D beam file. Relative paths are interpreted relative to the `constant/` directory.

## Tabulated beam file format

The tabulated beam file is a headerless ASCII file. It must not contain comment lines.

The required format is:

```text
nx ny
x0 y0
dx dy
f00 f10 f20 ... f(nx-1,0)
f01 f11 f21 ... f(nx-1,1)
...
f0(ny-1) ... f(nx-1,ny-1)
```

where:

| Entry | Meaning |
|---|---|
| `nx`, `ny` | Number of grid points in the x and y directions |
| `x0`, `y0` | Minimum x and y coordinates of the table, in meters |
| `dx`, `dy` | Uniform grid spacing in x and y, in meters |
| `f` | Relative planar intensity values |

The intensity values are stored in row-major order with `i` varying fastest:

```text
f_[i + nx*j]
```

The table is interpreted as nodal data. Bilinear interpolation is used between nodes.

The valid interpolation domain is:

```text
x0 <= x <= x0 + (nx - 1)*dx
y0 <= y <= y0 + (ny - 1)*dy
```

Outside this domain, the heat source returns zero.

## Example table

The included file:

```text
constant/beamProfile.txt
```

contains a 5 micron-resolution tabulated profile made from four overlapping Gaussian rings. It uses:

```text
nx ny: 101 101
x0 y0: -250 um, -250 um
dx dy: 5 um, 5 um
domain: -250 um to +250 um in x and y
```

The four rings are centered at:

```text
(-45 um, -45 um)
( 45 um, -45 um)
(-45 um,  45 um)
( 45 um,  45 um)
```

Each ring has a radius of 88 um and a Gaussian width of 21 um. The center-to-center spacing between neighboring rings is 90 um, so the rings overlap strongly and form a connected four-lobed beam profile.

The table is normalized so that the planar integral is approximately 1.0. This normalization is convenient, but the scan-path laser power still controls the total applied power.

The tabulated file defines only the relative 2D beam shape. The laser power should still be set in `constant/scanPath`, exactly as with the other heat source models.
