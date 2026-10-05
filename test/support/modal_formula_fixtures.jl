# A user-owned modal formula with fixed modal-to-phase maps. Items that need it take this
# module, so the items that use only `FormulaFixtures` do not reach ModalAnalysis.
@testmodule ModalFormulaFixtures begin
    using LineCableModels

    struct FixedModalMaps{P, O} <: LineCableModels.AbstractFormulation
        parameters::P
        options::O
    end
    FixedModalMaps(voltage::AbstractArray{<:Number, 3},
        current::AbstractArray{<:Number, 3}) = FixedModalMaps(
        (voltage = voltage, current = current), FormulationOptions())
    LineCableModels.formulation_options(selected::FixedModalMaps) = selected.options
    LineCableModels.formula_id(::FixedModalMaps) = :FixedModalMaps
    LineCableModels.formula_id(::Type{<:FixedModalMaps}) = :FixedModalMaps
    LineCableModels.description(::FixedModalMaps; compact::Bool = false) = "FixedModalMaps"
    LineCableModels.description(::Type{<:FixedModalMaps}; compact::Bool = false) =
        "FixedModalMaps"
    Base.NamedTuple(selected::FixedModalMaps) = (identifier = :FixedModalMaps,
        parameters = selected.parameters, options = selected.options.data)
    LineCableModels.description(
        ::Type{<:FixedModalMaps}, ::Val{:voltage}, value::AbstractArray;
        compact::Bool = false) = "voltage map="*sprint(show, value)
    LineCableModels.description(
        ::Type{<:FixedModalMaps}, ::Val{:current}, value::AbstractArray;
        compact::Bool = false) = "current map="*sprint(show, value)
    LineCableModels.Commons.initialize_buffers(
        ::FixedModalMaps, ::Type, input, plan, buffers) = buffers
    function LineCableModels.ModalAnalysis.decompose!(selected::FixedModalMaps,
            workspace, parameters, options)
        copyto!(workspace.Tv, selected.parameters.voltage)
        copyto!(workspace.Ti, selected.parameters.current)
        for frequency in axes(workspace.roots, 2), mode in axes(workspace.roots, 1)

            workspace.roots[mode, frequency] = sqrt((workspace.input.Y[
                :, :, frequency] * workspace.input.Z[:, :, frequency])[
                mode, mode])*workspace.input.root_scale
        end
        return workspace
    end
end
