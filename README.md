# TropicalHomotopies.jl
This project implements a general framework for computing solving polynomial systems numerically, 
based on homotopies constructed from certain tropical stable intersections.

It uses functionality from Oscar.jl and HomotopyContinuation.jl.

The underlying ideas are to be explained in a forthcoming paper with Paul Helminck, Oskar Henriksson, and Yue Ren.

> [!WARNING]  
> This is an early version of the code, made public for demonstration purposes. 
> Feel free to use it (at own risk). Comments and bug reports are welcome. 

## Installation
The easiest way to get started using the contents of this is to manually download all files, and then run

```
include("PATH/TO/src/main.jl")
```
where `PATH/TO/src/main.jl` is replaced by the actual path to the file `src/main.jl`.
