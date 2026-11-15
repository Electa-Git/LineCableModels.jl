function images(interior, f, g)
    for (pass, image) in enumerate((interior.impedance, interior.potential))
        z = pass == 1 ? f(image) : g(image)
        if pass == 1
            fill!(image, z)
        end
    end
end
ranged(f, g) = [pass == 1 ? f() : g() for pass in 1:2]
function flagged(A, B)
    for (magnetic, image) in ((true, A), (false, B))
        magnetic && fill!(image, 0)
        x = !magnetic ? 1 : 2
    end
end
function recorded(A, B)
    for field in ((image = A, magnetic = true), (image = B, magnetic = false))
        field.magnetic ? fill!(field.image, 0) : fill!(field.image, 1)
        field[2] || continue
    end
end
function destructured(A, B)
    for (; image, magnetic) in ((image = A, magnetic = true), (image = B, magnetic = false))
        if magnetic
            fill!(image, 0)
        end
    end
end
function identified(interior)
    for image in (interior.impedance, interior.potential)
        image === interior.impedance || continue
        image == interior.potential && fill!(image, 0)
        image in (interior.impedance,) && fill!(image, 1)
        image !== first(interior) && fill!(image, 2)
    end
end
function assigned(A, B)
    fields = (A, B)
    for (pass, image) in enumerate(fields)
        pass > 1 && fill!(image, 0)
    end
end
const FIELDS = (sin, cos)
positioned(x) = [k == 2 ? f(x) : 2f(x) for (k, f) in pairs(FIELDS)]
zipped(images) = [flag ? image : -image for (flag, image) in zip((true, false), images)]
predicated(xs) = [isodd(pass) ? x : -x for pass in 1:2, x in xs]
chosen(a, b) = [ifelse(pass == 1, a, b) for pass in 1:2]
shifted(a, b) = [pass - 1 == 0 ? a : b for pass in eachindex((a, b))]
weighted(a, b) = [(pass == 1) * a + (pass != 1) * b for pass in 1:2]
sieved(f, A, B) = sum(f(image) for (pass, image) in enumerate((A, B)) if pass == 1)
function crossed(M)
    for i in 1:2, j in 1:2

        i == j && continue
        M[i, j] = 0
    end
end
reversed(A, B) = [pass < 2 ? A : B for (pass, image) in enumerate(reverse([A, B]))]
function combined(interior, inner)
    for (pass, image) in enumerate((interior.impedance, interior.potential))
        if isempty(inner) && pass == 1
            fill!(image, 0)
        elseif pass == 2 || isempty(image)
            fill!(image, 1)
        end
    end
end
tupled(a, b) = ntuple(pass -> pass == 1 ? a : b, 2)
function eachone(A, B)
    map((true, false), (A, B)) do magnetic, image
        magnetic ? fill!(image, 0) : fill!(image, 1)
    end
end
summed(f, g, x) = sum(k -> k == 1 ? f(x) : g(x), 1:2)
reflected(image, magnetic) = magnetic ? fill!(image, 0) : fill!(image, 1)
both(A, B) = (reflected(A, true), reflected(B, false))
optioned(x; mode = 1) = mode == 1 ? x : -x
called(x) = optioned(x; mode = 2) + optioned(x, mode = 1)
function ordered(x, pass::Int = 1)
    if pass == 1
        return x
    end
    return -x
end
twice(x) = ordered(x) + ordered(x, 2)
function keyworded(f, A, B)
    for (pass, image) in enumerate((A, B))
        f(image; pass = pass, image)
        pass == 2 && fill!(image, 0)
    end
end
