module ImportsProbeClean
public Owner, Consumer
module Owner
export exported_value
public public_value
exported_value() = 1
public_value() = 2
end
module Consumer
import ..Owner
using ..Owner: exported_value
explicit() = exported_value()
qualified() = Owner.public_value()
end
end
