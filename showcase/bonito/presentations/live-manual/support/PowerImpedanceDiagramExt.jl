module PowerImpedanceDiagramExt

# Showcase-local extraction of the direct GraphMakie projection and renderer
# proposed in PowerImpedance.jl MR 36 at commit
# 71183fe1251cd178bc6c8704594d197d4d988414. The PlotBuilder integration is
# intentionally excluded. Remove this module when the upstream extension is
# available from the public PowerImpedance source.

import GraphMakie
import Graphs
import Makie
import NetworkLayout
import PowerImpedance
import PowerImpedance: PowerFlowResult

const NetworkBuilder = PowerImpedance.NetworkBuilder

include("PowerImpedanceDiagramExt/projection.jl")
include("PowerImpedanceDiagramExt/render.jl")

function diagram(
        network::NetworkBuilder.NetworkState;
        layout = NetworkLayout.Stress(; seed = 1),
        positions = nothing,
        interactive::Bool = true,
        style::NamedTuple = (;)
)
    return _render_diagram(
        project_diagram(network);
        layout,
        positions,
        interactive,
        style
    )
end

function diagram(
        network::NetworkBuilder.NetworkState,
        powerflow::PowerFlowResult;
        layout = NetworkLayout.Stress(; seed = 1),
        positions = nothing,
        interactive::Bool = true,
        style::NamedTuple = (;)
)
    return _render_diagram(
        project_diagram(network, powerflow);
        layout,
        positions,
        interactive,
        style
    )
end

export diagram, update_positions!

end
