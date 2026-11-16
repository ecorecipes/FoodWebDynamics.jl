# FoodWebDynamics.jl

```@docs
FoodWebDynamics
```

The model follows SI §8 of Domínguez-García et al. (2019). The constants the SI leaves open,
and the places where it contradicts its own dataset, are settled against their Dataset S1:
- the metabolic exponent −0.25;
- no metabolic loss for plants;
- the handling-time constant h₀ = 0.4;
- the survival rule `N ≥ floor(0.9 iniN)`;
- diet weights kept after extinctions.

See the docstring of [`BioenergeticParameters`](@ref).
