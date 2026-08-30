routes(::Val{:Deri}) = (;)
assumptions(::Val{:Deri}) = (;)
propagation(::Val{:Deri}) = Val(:backend)
description(::Formula{:Deri}) = "Deri-Semlyen"

:Deri
