module TagsProbePlanted
module Laws
struct Law{ID} end
law(::Val{ID}) where {ID} = Law{ID}()
pick(kind::Symbol) = law(Val(kind))
choose(kind::Symbol; scale = 1) = (law(Val(kind)), scale)
inner(law::Law) = (law, Val(:inner))
struct Knob{ID} end
(knob::Knob{:dial})(::Val{:J}, x) = x
end
module Uses
import ..Laws: pick
pick(kind::String) = pick(Val(Symbol(kind)))
end
end
