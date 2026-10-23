module Alpha
import ..Commons
public area, radius, half
area(x) = x^2
_half(x) = x / 2
half(x) = x / 2
radius(x) = sqrt(x) / 2
depth(x) = x - 1
function outer(x)
    inner(y) = y + 1
    return inner(x)
end
function scaled(x)
    double = y -> 2y
    local shifted::Function = function (y)
        y + 1
    end
    return shifted(double(x))
end
function summed(xs)
    return sum(xs) do x
        function shift(y)
            return y + 1
        end
        return shift(x)
    end
end
struct _Absent end
const _ABSENT = _Absent()
struct Shapes end
width(::Shapes, x) = x
macro twice(ex)
    return :(2 * $(esc(ex)))
end
alias(x; scale = 1) = Commons.shared(x; scale = scale)
end
