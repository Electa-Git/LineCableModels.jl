module TagsProbePlanted
module Laws
struct Law{ID} end
law(::Val{ID}) where {ID} = Law{ID}()
pick(kind::Symbol) = law(Val(kind))
choose(kind::Symbol; scale = 1) = (law(Val(kind)), scale)
inner(law::Law) = (law, Val(:inner))
end
end
