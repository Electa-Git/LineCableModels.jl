routes(::Val{:DeriSemlyen1981}) = (;)
assumptions(::Val{:DeriSemlyen1981}) = (;)
propagation(::Val{:DeriSemlyen1981}) = Val(:backend)
"""
$(TYPEDSIGNATURES)

**Identification.** Complex ground-return-plane approximation for overhead
conductors. Deri, Tevan, Semlyen, and Castanheira supplied the analytical
justification and error analysis for the complex-depth equations proposed by
Dubanton and published by Gary.

**Expression.** For the complex earth penetration depth
``p=(j\\omega\\mu_0\\sigma_g)^{-1/2}``,

```math
Z_{e,ii}=\\frac{j\\omega\\mu_0}{2\\pi}
\\ln\\frac{2(h_i+p)}{r_i},
```

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}
\\ln\\frac{\\sqrt{(h_i+h_j+2p)^2+y_{ij}^2}}
{\\sqrt{(h_i-h_j)^2+y_{ij}^2}}.
```

**Reference.** A. Deri, G. Tevan, A. Semlyen, and A. Castanheira, “The
Complex Ground Return Plane: A Simplified Model for Homogeneous and Multi-Layer
Earth Return,” *IEEE Transactions on Power Apparatus and Systems*, PAS-100(8),
3686–3693, 1981. DOI: 10.1109/TPAS.1981.317011.
"""
description(::Formula{:DeriSemlyen1981}) = "Deri-Semlyen"

:DeriSemlyen1981
