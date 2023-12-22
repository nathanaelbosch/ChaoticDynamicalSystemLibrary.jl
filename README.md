# ChaoticDynamicalSystemLibrary.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://nathanaelbosch.github.io/ChaoticDynamicalSystemLibrary.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://nathanaelbosch.github.io/ChaoticDynamicalSystemLibrary.jl/dev/)
[![Build Status](https://github.com/nathanaelbosch/ChaoticDynamicalSystemLibrary.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/nathanaelbosch/ChaoticDynamicalSystemLibrary.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/nathanaelbosch/ChaoticDynamicalSystemLibrary.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/nathanaelbosch/ChaoticDynamicalSystemLibrary.jl)

ChaoticDynamicalSystemLibrary.jl is a collection of chaotic dynamical systems for the Julia [SciML](https://sciml.ai/) ecosystem.

The package is a partial port of the [`dysts`](https://github.com/williamgilpin/dysts) Python package developed by [William Gilpin](https://github.com/williamgilpin).
It implements many (but not all) of the same systems, but does not provide any functionality to simulate them.
In contrast to [`dysts`](https://github.com/williamgilpin/dysts), the main focus of this package is not to be a benchmark for general time-series ML models, but only to provide a collection of ODEs.
Their simulation is left to other packages, such as [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/).

---

**[**William Gilpin**](https://github.com/williamgilpin) deserves most of the credit for this package.
He is the original author of the [`dysts`](https://github.com/williamgilpin/dysts), and without [`dysts`](https://github.com/williamgilpin/dysts) this package would not exist.**

---



## Installation

```julia
import Pkg
Pkg.add("ChaoticDynamicalSystemLibrary")
```

## Usage

```julia
using ChaoticDynamicalSystemLibrary, DifferentialEquations, Plots

prob = ChaoticDynamicalSystemLibrary.Lorenz()
sol = solve(prob, Tsit5(), tspan=(0, 100))
plot(sol, idxs=(1, 2, 3))
```
![Lorenz Simulation](./docs/src/readmeplot.svg?raw=true "Lorenz Simulation")

For a full list of the available systems, see the [documentation](https://nathanaelbosch.github.io/ChaoticDynamicalSystemLibrary.jl/stable/).


## Acknowledgements

- [William Gilpin](https://github.com/williamgilpin) / [`dysts`](https://github.com/williamgilpin/dysts):
  The foundation that this package is based on.
- [Jürgen Meier](http://www.3d-meier.de/tut19/Seite1.html) and
  [J. C. Sprott](http://sprott.physics.wisc.edu/sprott.htm):
  `dysts` contains systems from both collections, and therefore so does this package.
