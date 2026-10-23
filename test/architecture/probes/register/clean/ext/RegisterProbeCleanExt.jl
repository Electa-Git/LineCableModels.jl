module RegisterProbeCleanExt
import ..RegisterProbeClean: Alpha
public draw
Alpha.area(x::Int) = x * x
draw(item) = Alpha.area(item) + 1
end
