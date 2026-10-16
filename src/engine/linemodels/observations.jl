Units.quantity(::typeof(poles)) = Units.Quantity{:pole}()
Units.quantity(::typeof(residues)) = Units.Quantity{:residue}()
Units.quantity(::typeof(delays)) = Units.Quantity{:delay}()
Units.quantity(::typeof(fitting_error)) = Units.Quantity{:fitting_error}()

"""Select the RMS relative fitting-error matrix for one quantity."""
fitting_error(result::LineModel, quantity::Symbol) = getproperty(result.fitting_error, quantity).rms

observe(result::LineModel, ::typeof(poles)) = poles(result)
observe(result::LineModel, ::typeof(delays)) = delays(result)
observe(result::LineModel, ::typeof(fitting_error), quantity::Symbol) = fitting_error(result, quantity)

observables(::Type{<:LineModel}) = (poles, delays, fitting_error)
basis(::LineModel) = :total

function Commons.observation_requests(result::LineModel, requests::Tuple; complete_pairs::Bool=false)
    selected = isempty(requests) ? observables(typeof(result)) : requests
    retained = Tuple(Iterators.flatten(request === fitting_error ?
        ((fitting_error, :admittance), (fitting_error, :propagation)) :
        (request,) for request in selected))
    Commons.validate_observables(result, retained, ())
    allunique(retained) || throw(ArgumentError("observation requests must be distinct"))
    return (retained=retained, displayed=selected)
end

function Commons.observation_quantity(result::LineModel, request::Tuple{typeof(fitting_error), Symbol};
        unit=nothing, clip=true, atol=nothing, frequencies=nothing)
    resolved = Commons.observation_request(result, request)
    raw = observe(result, request...)
    q = resolved.quantity
    native = Units.native_unit(q, basis(result))
    target = unit === nothing ? Units.display_unit(q, basis(result)) : unit
    return (request, quantity=q, family=last(request), statistic=:rms,
        values=Commons.detach(raw, Units.scale_factor(native, target)), unit=target,
        basis=basis(result), coordinates=(kind=:array, indices=resolved.indices,
            extent=size(raw), labels=copy(result.coordinates)), thresholds=nothing,
        available=nothing, engineering_zero=nothing, clipped=false, missing_reason=nothing)
end
