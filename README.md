# TropicalHomotopies.jl
This project implements a general framework for computing solving polynomial systems numerically, 
based on homotopies constructed from certain tropical stable intersections.

It uses functionality from Oscar.jl and HomotopyContinuation.jl.

The underlying ideas are to be explained in a forthcoming paper with Paul Helminck, Oskar Henriksson, Yue Ren, and Benjamin Schröter.

> [!WARNING]  
> This is an early version of the code, made public for demonstration purposes. 
> Feel free to use it (at own risk). Comments and bug reports are welcome. 

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
### Principal decomposition
The following decomposition traces the mixed volume number of paths:
```
decomposition = [[f] for f in F]
tropical_solve(decomposition)
```
The underlying stable intersection calculation shows that the mixed volume is 4:
```
tropical_root_bound(decomposition)
```
### Vertical decomposition
A sharper bound is obtained by using the linear--binomial ("vertical") decomposition, as follows:
```
decomposition = linear_and_binomial_part(F)
tropical_root_bound(decomposition)
```
To use the associated generalized polyhedral homotopies to solve the system, we run
```
tropical_solve(decomposition)
```
There is also a wrapper for these functions, that also includes a check for transversality of the matroid of the linear part, which allows for a faster solving via traditional polyhedral homotopies:
```
vertical_solve(F)
```
