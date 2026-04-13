## The short version

1.  EMC2 only implements the fixed-boundary versions of these models and no `s_t`.
2.  EMC2 parameterizes drift by magnitude and angles, whereas CRDDM uses Cartesian drift vectors.
3.  EMC2 does not assume `sigma = 1`
4.  EMC2 simulates boundary hits with interpolation.
5.  In the spherical and hyperspherical models, EMC2 treats the response likelihood as a density in the observed angular coordinates.
6.  Projected versions are made as marginals of their non-projected.

## Family-wide differences

### 1. Parameters

EMC2 does not have any of the threshold functions nor s_t.

### 2. Drift parameterization

CRDDM works with Cartesian drift vectors:

-   `Circular`: 2D drift vector
-   `Spherical`: 3D drift vector
-   `HyperSpherical`: 4D drift vector

EMC2 instead parameterizes drift as:

-   a drift magnitude `v`
-   one or more angular coordinates (polar and azimuth) `theta`, `theta1`, `theta2`, `theta3`

Theta and theta1 is not used very consistently. Maybe make it theta? Maybe theta1?

### 3. `sigma` handling

In EMC2, the fixed-boundary likelihoods consistently use the scaled time

``` text
t sigma^2 / a^2
```

and carry `sigma^2` through the drift correction and `sv` integration terms.

CRDDM uses `sigma` in the first-passage-time expressions, but several other fixed-boundary terms assume `sigma = 1`.

### 4. Simulation of the boundary hit

CRDDM stops when the process first lands outside the boundary and uses the end of that step as the observed hit.

In EMC2: after detecting that the final step crossed the boundary, it interpolates along that last segment to estimate the actual hit point on the boundary. That changes both the RT and the response angle. I found this matters for SBC. However, I also increased number of simulation steps, which makes it a bit slow to simulate at the moment. It might well be that a smaller number is enough

## Model by model

## CDM

No differences outside of the above.

## SDM

CRDDM writes the angular part as a density on the sphere itself. EMC2 writes it as a density in the angular coordinates that are actually observed and returned by the simulator, namely `R` and `R2`.

On a sphere,

``` text
dA = sin(R) dR dR2
```

so a density with respect to surface area is not the same thing as a density with respect to raw angle coordinates. To convert from one to the other, a factor of `sin(R)` is needed.

That factor is present in EMC2 and absent in CRDDM.

So the current EMC2 `SDM` is the normalized density in the observed variables `R` and `R2`. CRDDM’s `Spherical` likelihood is closer to a density with respect to surface measure.

## HSDM

Same here, the response lives on the 3-sphere, and the hyperspherical coordinate Jacobian is

``` text
sin(R)^2 sin(R2)
```

So if the observed responses are `R`, `R2`, and `R3`, I think the coordinate-density needs that factor.

There is also a second difference in `HSDM`: the normalization constant. EMC2 uses what I think is the correct surface-area constant for the 3-sphere,

``` text
2*pi^2
```

whereas CRDDM uses `2*pi`.

## PSDM

In EMC2, `PSDM` is built as the marginal of the full `SDM`. The projected kernel has the modified-Bessel form, and the normalization follows from integrating the parent model over the hidden angle.

There are really two differences here:

1.  The projected formula in EMC2 is internally tied to the full `SDM` by construction.
2.  EMC2 then interprets the result as a density in the observed angle `R`, which adds the same coordinate Jacobian as above.

On the sphere, that means the observable density in `R` carries a factor of `sin(R)`.

## PHSDM

As in `PSDM`, EMC2’s implementation behaves like a marginal of the full parent model. The projected kernel has the modified-Bessel structure, and the normalization is consistent with integrating over the hidden final angle.

EMC2’s `PHSDM` also includes the hyperspherical coordinate Jacobian for the observed variables `R` and `R2`, namely

``` text
sin(R)^2 sin(R2)
```

In EMC2, the projected models are best understood as the parent models with one angle integrated out.
