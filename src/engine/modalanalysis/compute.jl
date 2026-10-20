function offdiagonal_ratio(matrix::AbstractMatrix)
    n = checksquare(matrix)
    n > 0 || throw(ArgumentError("modal matrix must be nonempty"))
    first_value = nominal(matrix[1, 1])
    matrix_norm_squared = zero(real(abs2(first_value)))
    offdiagonal_norm_squared = zero(matrix_norm_squared)
    @inbounds for column in 1:n, row in 1:n
        squared = abs2(nominal(matrix[row, column]))
        matrix_norm_squared += squared
        row == column || (offdiagonal_norm_squared += squared)
    end
    matrix_norm = sqrt(matrix_norm_squared)
    return sqrt(offdiagonal_norm_squared) /
           max(matrix_norm, eps(float(real(matrix_norm))))
end


# The modal backend normalizes the computation options of a modal-analysis formulation.
computation_options(::Type{<:ModalAnalysisFormulation}, record::ComputationOptions) =
    computation_options(LineCableModelsModal, record)

# A normalized record is returned as it is.
computation_options(::Type{LineCableModelsModal},
    record::ComputationOptions{<:NamedTuple{(:offdiagonal_tolerance, :on_result, :timing,
        :rotate, :verbosity), <:Tuple{Real, Any, Bool, Bool, Any}}}) = record

function computation_options(::Type{LineCableModelsModal}, record::ComputationOptions)::ComputationOptions
    options = record.data
    isempty(setdiff(keys(options), (:rotate, :offdiagonal_tolerance, :on_result, :timing, :verbosity))) ||
        throw(ArgumentError("unknown modal action options: $(Tuple(keys(options)))"))
    tolerance = get(options, :offdiagonal_tolerance, 1e-6)
    tolerance isa Real && isfinite(tolerance) && tolerance >= 0 ||
        throw(ArgumentError("offdiagonal_tolerance must be finite and nonnegative"))
    callback = get(options, :on_result, nothing)
    (callback === nothing || callback isa Function) ||
        throw(ArgumentError("on_result must be a callable function or nothing"))
    timing = get(options, :timing, false)
    timing isa Bool || throw(ArgumentError("timing must be Bool"))
    rotate = get(options, :rotate, true)
    rotate isa Bool || throw(ArgumentError("rotate must be Bool"))
    return ComputationOptions(; offdiagonal_tolerance=tolerance, on_result=callback,
        timing, rotate, verbosity=get(options, :verbosity, 0))
end

"""Store one scalar modal computation's results, common arrays, selected work, and diagnostics."""
struct ModalAnalysisWorkspace{P,I,V,R,Z,Y,B,D}
    source::P
    input::I
    Tv::V
    Ti::V
    roots::R
    Zm::Z
    Ym::Y
    buffers::B
    diagnostics::D
end

formula_parameters(selected::Formula) = selected.parameters
formula_parameters(::AbstractFormulation) = (;)

# Total coefficients scale the roots by the line length. Per unit length, the scale is `unit`.
_root_scale(::Units.PerUnitLength, ell0, unit) = unit
_root_scale(::Units.TotalLength, ell0, unit) = nominal(ell0)

function ModalAnalysisWorkspace(source::LineParameters, selected::AbstractFormulation)
    Zp = nominal(source.Z.values)
    Yp = nominal(source.Y.values)
    n, _, nf = size(Zp)
    ell0 = line_length(source)
    factor = _root_scale(basis(source), ell0, one(real(zero(eltype(Zp)))))
    input = (Z=Zp, Y=Yp, f=nominal(source.f),
        root_scale=factor, source_basis=basis(source))
    T = float(promote_type(eltype(input.Z), eltype(input.Y)))
    Tv = Array{T,3}(undef, n, n, nf)
    Ti = similar(Tv)
    roots = Array{T,2}(undef, n, nf)
    S = promote_type(eltype(source.Z.values), eltype(Tv))
    Zm = Array{S,3}(undef, n, n, nf)
    Ym = similar(Zm)
    plan = (; n, nf)
    common = (Zslice=Matrix{T}(undef,n,n),Yslice=Matrix{T}(undef,n,n),
        coordinate_product=Matrix{S}(undef,n,n),
        coordinate_factor=Matrix{T}(undef,n,n),
        admittance_impedance_product=Matrix{T}(undef,n,n),
        previous_eigenvalues=Vector{T}(undef,n),eigenvalues=Vector{T}(undef,n),
        previous_eigenvectors=Matrix{T}(undef,n,n),eigenvectors=Matrix{T}(undef,n,n),
        propagation_eigenvalues=Vector{T}(undef,n),voltage_vector=Vector{T}(undef,n))
    buffers = initialize_buffers(selected, T, input, plan, common)
    R=typeof(real(zero(T)))
    residual=Matrix{Union{Nothing,R}}(undef,n,nf)
    iteration_counts=Matrix{Union{Nothing,Int}}(undef,n,nf)
    convergence=Matrix{Union{Nothing,Bool}}(undef,n,nf)
    fill!(residual,nothing)
    fill!(iteration_counts,nothing)
    fill!(convergence,nothing)
    diagnostics = (fallback_frequencies=Int[], missed_frequencies=Int[],
        z_coupling=Vector{typeof(real(zero(T)))}(undef, nf),
        y_coupling=Vector{typeof(real(zero(T)))}(undef, nf),
        eigen_residual=residual,iterations=iteration_counts,converged=convergence)
    return ModalAnalysisWorkspace(source, input, Tv, Ti, roots, Zm, Ym, buffers, diagnostics)
end

function decompose! end

function _coordinate_algebra!(workspace::ModalAnalysisWorkspace, execution)
    source = workspace.source
    nfreq = size(workspace.Zm, 3)
    tolerance = execution.data.offdiagonal_tolerance
    product=workspace.buffers.coordinate_product
    factor=workspace.buffers.coordinate_factor
    for frequency in 1:nfreq
        Tvf = @view workspace.Tv[:, :, frequency]
        Tif = @view workspace.Ti[:, :, frequency]
        Zp = @view source.Z.values[:, :, frequency]
        Yp = @view source.Y.values[:, :, frequency]
        mul!(product,Zp,Tif)
        copyto!(factor,Tvf)
        ldiv!(lu!(factor),product)
        copyto!(@view(workspace.Zm[:,:,frequency]),product)
        mul!(product,Yp,Tvf)
        copyto!(factor,Tif)
        ldiv!(lu!(factor),product)
        copyto!(@view(workspace.Ym[:,:,frequency]),product)
        zr = offdiagonal_ratio(@view(workspace.Zm[:, :, frequency]))
        yr = offdiagonal_ratio(@view(workspace.Ym[:, :, frequency]))
        workspace.diagnostics.z_coupling[frequency] = zr
        workspace.diagnostics.y_coupling[frequency] = yr
        if zr > tolerance || yr > tolerance
            @warn "modal coefficients exceed the requested off-diagonal tolerance" frequency zr yr tolerance
        end
    end
    return workspace
end

# Differentiate the diagonal product at the retained nominal branch. The maps
# stay nominal. Uncertainty in the transformed coefficients remains connected
# to its original scalar sources.
function _dependent_roots(workspace::ModalAnalysisWorkspace)
    nominal_roots=workspace.roots
    eltype(workspace.Zm)===eltype(nominal_roots) && return copy(nominal_roots)
    S=promote_type(eltype(workspace.Zm),eltype(nominal_roots))
    roots=Array{S}(undef,size(nominal_roots))
    for frequency in axes(roots,2), mode in axes(roots,1)
        root=nominal_roots[mode,frequency]
        iszero(root) && throw(DomainError(root,
            "first-order modal root is undefined at a zero nominal root"))
        product=workspace.Zm[mode,mode,frequency]*workspace.Ym[mode,mode,frequency]
        roots[mode,frequency]=root+(product-nominal(product))/(2root)
    end
    return roots
end

"""
$(TYPEDSIGNATURES)

Read the modal `formulation` that a point applied to the result whose record is `source`, as
data. The formulations and selections are those of the source. The descriptions list the
source's selections, then the modal transformation, for every quantity.
"""
function Engine.completed_formulation(formulation::ModalAnalysisFormulation,
        source::ComputationDetails)
    upstream=Engine.completed_formulation(source)
    modal_descriptions=filter(field -> !isempty(field.meaning),
        Engine.completed_formulation(formulation).formulation_fields.all)
    descriptions=vcat(get(upstream.formulation_fields,:all,NamedTuple[]),modal_descriptions)
    return (formulations=upstream.formulations,selections=upstream.selections,
        formulation_fields=merge(upstream.formulation_fields,
            (all=descriptions,Z=copy(descriptions),Y=copy(descriptions))))
end

"""
$(TYPEDSIGNATURES)

Compute the modal analysis of one point: the modal decomposition of the phase-domain line
parameters of the problem.
"""
function compute(
        resolved::Resolved{<:ModalAnalysisProblem{<:LineParameters{T,U,PhaseDomain,Basis}}},
        formulation::ModalAnalysisFormulation;
        options::Union{NamedTuple,ComputationOptions}=ComputationOptions()) where {T,U,Basis}
    execution = computation_options(LineCableModelsModal, ComputationOptions(options))
    parameters=resolved.problem.parameters
    selected = formulation.formula
    workspace = ModalAnalysisWorkspace(parameters, selected)
    decompose!(selected, workspace, formula_parameters(selected),
        formulation_options(selected))
    all(isfinite, workspace.Tv) || throw(DomainError(workspace.Tv, "Tv must be finite"))
    all(isfinite, workspace.Ti) || throw(DomainError(workspace.Ti, "Ti must be finite"))
    orient_modes!(workspace.Tv, workspace.Ti, execution.data.rotate)
    maps = ModalOperators(copy(workspace.Tv), copy(workspace.Ti))
    _coordinate_algebra!(workspace, execution)
    roots = _dependent_roots(workspace)
    all(isfinite, roots) || throw(DomainError(roots, "modal roots must be finite"))
    modal = ModalDomain(maps, roots)
    diagnostics = map(copy,workspace.diagnostics)
    # The record keeps the source's entries, except its timing, then adds the modal entries,
    # among them the source's record as the choice of the modal formulation.
    source = parameters.details
    kept = filter(∉((:timing, :comparison_unsupported, :source_gridpoint, :modal, :formulation,
        :choices)), keys(source.data))
    retained = ComputationDetails(merge(NamedTuple{kept}(source.data),
        NamedTuple{(:source_gridpoint, :modal, :formulation, :choices),
            Tuple{Any, NamedTuple, ModalAnalysisFormulation, ComputationDetails}}((
            get(source.data, :gridpoint, nothing),
            (identifier=formula_id(selected), requested=NamedTuple(formulation.definition),
                effective=NamedTuple(selected), rotate=execution.data.rotate,
                diagnostics=diagnostics),
            formulation, source))))
    return LineParameters(modal,
        SeriesImpedance{eltype(workspace.Zm),Basis}(copy(workspace.Zm)),
        ShuntAdmittance{eltype(workspace.Ym),Basis}(copy(workspace.Ym)),
        parameters.f, parameters.terminals, parameters.line_length, retained)
end

"""Restore phase coefficients using the retained modal-to-phase bases."""
function transform(::Type{PhaseDomain}, parameters::LineParameters{T,U,D,Basis}) where {T,U,D<:ModalDomain,Basis}
    maps = parameters.domain.operators
    dimensions = size(parameters.Z.values)
    S = promote_type(T, eltype(maps.Tv), eltype(maps.Ti))
    impedance = Array{S,3}(undef, dimensions)
    admittance = similar(impedance)
    for frequency in axes(impedance, 3)
        Tvf = @view maps.Tv[:, :, frequency]
        Tif = @view maps.Ti[:, :, frequency]
        @views impedance[:, :, frequency] .= (Tvf * parameters.Z.values[:, :, frequency]) / Tif
        @views admittance[:, :, frequency] .= (Tif * parameters.Y.values[:, :, frequency]) / Tvf
    end
    return LineParameters(PhaseDomain(),
        SeriesImpedance{S,Basis}(impedance),
        ShuntAdmittance{S,Basis}(admittance), parameters.f,
        parameters.terminals, parameters.line_length, parameters.details)
end

description(::ModalAnalysisProblem; compact::Bool = false) = "Modal analysis"

function computation_details(::Type{<:ModalAnalysisFormulation}, result::LineParameters)::ComputationDetails
    return details(result)
end
