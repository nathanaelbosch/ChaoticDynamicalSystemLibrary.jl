```@meta
CurrentModule = ChaoticDynamicalSystemLibrary
```

# ChaoticDynamicalSystemLibrary.jl

ChaoticDynamicalSystemLibrary.jl is a collection of chaotic dynamical systems for the Julia [SciML](https://sciml.ai/) ecosystem.

The package is a partial port of the [`dysts`](https://github.com/williamgilpin/dysts) Python package developed by [William Gilpin](https://github.com/williamgilpin).
It implements many (but not all) of the same systems, but does not provide any functionality to simulate them.
In contrast to [`dysts`](https://github.com/williamgilpin/dysts), the main focus of this package is not to be a benchmark for general time-series ML models, but only to provide a collection of ODEs.
Their simulation is left to other packages, such as [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/).

---

**[William Gilpin](https://github.com/williamgilpin) deserves most of the credit for this package.
He is the original author of the [`dysts`](https://github.com/williamgilpin/dysts), and without [`dysts`](https://github.com/williamgilpin/dysts) this package would not exist.**

---

## Full list of all dynamical systems:

```@autodocs
Modules = [ChaoticDynamicalSystemLibrary]
Public = false
Pages   = ["chaotic_attractors.jl"]
```
