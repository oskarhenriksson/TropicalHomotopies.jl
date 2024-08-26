# TropicalHomotopies.jl
This project implements a general framework for computing solving polynomial systems numerically, 
based on homotopies constructed from certain tropical stable intersections.

It uses functionality from Oscar.jl and HomotopyContinuation.jl.

The underlying ideas are to be explained in a forthcoming paper with Paul Helminck, Oskar Henriksson, and Yue Ren.

> [!WARNING]  
> This is an early version of the code, made public for demonstration purposes. 
> Feel free to use it (at own risk). Comments and bug reports are welcome. 

## Structure of the repository
The repository consists of three main parts:
* The `examples` directory contains explicit computations used for the examples given throughout the paper.
* The `src` directory contains an implementation of the main algorithm from the paper in the case of vertically parametrized systems.
* The `case_studies` contains code for the examples from chemical reaction network theory and graph rigidity theory discussed at the end of the paper, that uses functionality from `src`.

## Installation
The easiest way to get started using the contents of this is to manually download all files, and then run

```
include("PATH/TO/src/main.jl")
```
where `PATH/TO/src/main.jl` is replaced by the actual path to the file `src/main.jl`.
