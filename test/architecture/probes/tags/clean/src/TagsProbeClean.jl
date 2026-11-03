module TagsProbeClean
module Laws
struct Law{ID} end
Law(::Val{ID}) where {ID} = Law{ID}()
Law(kind::Symbol; strict = true) = Law(Val(kind))
law(::Val{ID}) where {ID} = Law{ID}()
law(kind::Symbol) = law(Val(kind))
part(law::Law, ::Val{:inner}) = law
inner(law::Law) = part(law, Val(:inner))
struct Formula{ID} end
(formula::Formula{:linear})(::Val{:J}, x) = x
(formula::Formula{:linear})(::Val{:Z}, ::Val{:self}, i, x) = x
end
module Uses
import ..Laws: law
law(kind::String) = law(Val(Symbol(kind)))
end
end
