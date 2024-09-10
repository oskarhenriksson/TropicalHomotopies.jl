# TropicalHomotopies.jl
This project is a proof of concept a general framework for solving polynomial systems numerically, based on homotopies constructed using tropical geometry. It uses functionality from [OSCAR](https://github.com/oscar-system/Oscar.jl) and [HomotopyContinuation.j](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl).

The underlying ideas are to be explained in a forthcoming manuscript titled _A tropical method for solving parametrized polynomial systems_ by Paul Helminck, Oskar Henriksson, and Yue Ren.

## Structure of the repository
The repository consists of three main parts:
* The `examples` directory contains explicit computations used for the examples given throughout the paper.
* The `src` directory contains an implementation of the main algorithm from the paper in the case of vertically parametrized systems.
* The `case_studies` contains code for the examples from chemical reaction network theory and graph rigidity theory discussed at the end of the paper, that uses functionality from `src`.

> [!WARNING]  
> The path-tracking portion of the code is still experimental, and sometimes displays numerical instabilities.

## Examples

See the notebook `case_studies/wnt_pathway.ipynb`.

## Installation
The easiest way to get started using the code in `src` is to manually download all files, and then run

```julia
include("PATH/TO/src/main.jl")
```
where `PATH/TO/src/main.jl` is replaced by the actual path to the file `src/main.jl`.

> [!WARNING]  
> Parts of the code require the latest development version of OSCAR (09/09/2024). 
