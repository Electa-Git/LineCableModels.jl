module ImportsProbePlanted
public Exporter, Owner, Relay, Consumer
module Exporter
export implicit_value
implicit_value() = 1
end
module Owner
public relayed_value, public_value
relayed_value() = 2
public_value() = 3
private_value() = 4
private_import() = 5
end
module Relay
using ..Owner: relayed_value
relay() = relayed_value()
end
module Consumer
using ..Exporter
import ..Owner
import ..Relay
using ..Owner: public_value
using ..Owner: private_import
using ..Relay: relayed_value
implicit() = implicit_value()
relayed() = relayed_value()
imported() = private_import()
private() = Owner.private_value()
relay_access() = Relay.relayed_value()
self() = Consumer.implicit()
end
end
