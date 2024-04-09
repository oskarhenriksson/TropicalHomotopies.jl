# TropicalHomotopies.jl
This projects implement a general framework for computing solving polynomial systems numerically, 
based on homotopies constructed from certain tropical stable intersections.

It uses functionality from Oscar.jl and HomotopyContinuation.jl.

## Installation
The easiest way to get started using this is to download the project, and run

```
include(path_to_TropicalHomotopies.jl)
```

## Example
Consider the following system:
```
Qx, x = polynomial_ring(QQ,"x"=>1:3);
F = [5*x[1]^3*x[2] - 6*x[1]*x[2]^3 + x[1]*x[2], 
    5*x[1]^3*x[2] - 6*x[1]*x[2]^3 - x[1]*x[2]^2,
    x[1]+x[2]+x[3]-2]
```
We can solve it with classical polyhedral homotopies as follows:
```
systems = [[f] for f in F]
tropical_solve(systems)
```
The underlying stable intersection calculation shows that the mixed volume is 4:
```
tropical_stable_intersection_with_homotopy_data(systems)[1]
```
A sharper bound is obtained as follows:
```
systems = linear_and_binomial_part(F)
tropical_stable_intersection_with_homotopy_data(systems)[1]
```
and the corresponding homotopies are given by
```
tropical_solve(systems)
```
