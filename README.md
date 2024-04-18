# TropicalHomotopies.jl
This project implements a general framework for computing solving polynomial systems numerically, 
based on homotopies constructed from certain tropical stable intersections.

It uses functionality from Oscar.jl and HomotopyContinuation.jl.

## Installation
The easiest way to get started using the contents of this is to manually download all files, and then run

```
include(TropicalHomotopies.jl)
```
where `TropicalHomotopies.jl` is replaced by the path to the file `TropicalHomotopies.jl`.

## Examples
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
tropical_root_bound(systems)
```
A sharper bound is obtained as follows:
```
systems = linear_and_binomial_part(F)
tropical_root_bound(systems)
```
To use the associated generalized polyhedral homotopies to solve the system, we run
```
tropical_solve(systems)
```
