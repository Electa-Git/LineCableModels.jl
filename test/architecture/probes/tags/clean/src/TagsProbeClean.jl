module TagsProbeClean
module Laws
struct Law{ID} end
Law(::Val{ID}) where {ID} = Law{ID}()
Law(kind::Symbol; strict = true) = Law(Val(kind))
law(::Val{ID}) where {ID} = Law{ID}()
law(kind::Symbol) = law(Val(kind))
part(law::Law, ::Val{:inner}) = law
inner(law::Law) = part(law, Val(:inner))
end
end
