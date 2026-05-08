# TropicalHomotopies.jl
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://oskarhenriksson.github.io/TropicalHomotopies.jl/dev/)
[![CI](https://github.com/oskarhenriksson/TropicalHomotopies.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/oskarhenriksson/TropicalHomotopies.jl/actions/workflows/ci.yml)

This project is a proof of concept for a general framework for solving polynomial systems through homotopies constructed using tropical geometry, developed in the paper [A tropical method for solving parametrized polynomial systems](https://arxiv.org/abs/2409.13288) by Paul Helminck, Oskar Henriksson, and Yue Ren. The code uses functionality from [OSCAR](https://github.com/oscar-system/Oscar.jl) and [HomotopyContinuation.j](https://github.com/JuliaHomotopyContinuation/HomotopyContinuation.jl).

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
The easiest way to get started using the code is to run the following commands:

```julia
using Pkg
Pkg.add(url="https://github.com/oskarhenriksson/TropicalHomotopies.jl")
```