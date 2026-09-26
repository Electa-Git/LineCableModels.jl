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

function _check_operators(maps::ModalOperators, parameters::LineParameters)
    expected = size(parameters.Z.values)
    size(maps.Tv) == expected || throw(DimensionMismatch("Tv must have size $expected"))
    size(maps.Ti) == expected || throw(DimensionMismatch("Ti must have size $expected"))
    all(isfinite, maps.Tv) || throw(DomainError(maps.Tv, "Tv must be finite"))
    all(isfinite, maps.Ti) || throw(DomainError(maps.Ti, "Ti must be finite"))
    return maps
end

function computation_options(::Type{LineCableModelsModal}, record::ComputationOptions)::ComputationOptions
    options = record.data
    isempty(setdiff(keys(options), (:offdiagonal_tolerance, :on_result, :timing, :verbosity))) ||
        throw(ArgumentError("unknown modal action options: $(Tuple(keys(options)))"))
    tolerance = get(options, :offdiagonal_tolerance, 1e-6)
    tolerance isa Real && isfinite(tolerance) && tolerance >= 0 ||
        throw(ArgumentError("offdiagonal_tolerance must be finite and nonnegative"))
    callback = get(options, :on_result, nothing)
    (callback === nothing || callback isa Function) ||
        throw(ArgumentError("on_result must be a callable function or nothing"))
    timing = get(options, :timing, false)
    timing isa Bool || throw(ArgumentError("timing must be Bool"))
    return ComputationOptions(; offdiagonal_tolerance=tolerance, on_result=callback,
        timing, verbosity=get(options, :verbosity, 0))
end

"""Store one scalar modal calculation's results, common arrays, selected work, and diagnostics."""
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

allocation_selector(::Formula{ID}) where {ID} = Val(ID)
allocation_selector(selected::AbstractFormulation) = selected
formula_parameters(selected::Formula) = selected.parameters
formula_parameters(::AbstractFormulation) = (;)

function ModalAnalysisWorkspace(source::LineParameters, selected::AbstractFormulation)
    Zp = nominal(source.Z.values)
    Yp = nominal(source.Y.values)
    n, _, nf = size(Zp)
    ell0 = line_length(source)
    factor = one(real(zero(eltype(Zp))))
    if basis(source) === :total && ell0 !== nothing
        ell0 isa Real && !(ell0 isa Bool) && isfinite(nominal(ell0)) && nominal(ell0) > 0 ||
            throw(ArgumentError("total coefficients require a positive source normalization length"))
        factor = nominal(ell0)
    end
    input = (Z=Zp, Y=Yp, f=nominal(source.f),
        root_scale=factor, source_basis=basis(source))
    T = float(promote_type(eltype(input.Z), eltype(input.Y)))
    Tv = Array{T,3}(undef, n, n, nf)
    Ti = similar(Tv)
    roots = Array{T,2}(undef, n, nf)
    S = promote_type(eltype(source.Z.values), eltype(Tv))
    Zm = Array{S,3}(undef, n, n, nf)
    Ym = similar(Zm)
    invariants = (; n, nf)
    common = (Zslice=Matrix{T}(undef,n,n),Yslice=Matrix{T}(undef,n,n),
        coordinate_product=Matrix{S}(undef,n,n),
        coordinate_factor=Matrix{T}(undef,n,n),
        admittance_impedance_product=Matrix{T}(undef,n,n),
        previous_eigenvalues=Vector{T}(undef,n),eigenvalues=Vector{T}(undef,n),
        previous_eigenvectors=Matrix{T}(undef,n,n),eigenvectors=Matrix{T}(undef,n,n),
        propagation_eigenvalues=Vector{T}(undef,n),voltage_vector=Vector{T}(undef,n))
    buffers = initialize_buffers(allocation_selector(selected), T, input, invariants, common)
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
# stay nominal; uncertainty in the transformed coefficients remains connected
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

# Keep the completed record's outer key/type layout inferable even when an
# upstream gridpoint or formula description is intentionally type-erased.
@generated function _modal_detail_merge(record::NamedTuple{Names,Types},
        extra::NamedTuple{ExtraNames,ExtraTypes}) where {Names,Types,ExtraNames,ExtraTypes}
    retained=filter(name -> name ∉ ExtraNames &&
        name ∉ (:timing,:comparison_unsupported),Names)
    names=(retained...,ExtraNames...)
    types=(map(name -> fieldtype(Types,findfirst(==(name),Names)),retained)...,
        fieldtypes(ExtraTypes)...)
    compact=map(type -> type<:NamedTuple ? NamedTuple : type,types)
    output=NamedTuple{names,Tuple{compact...}}
    entries=[:(getproperty(record,$(QuoteNode(name)))) for name in retained]
    append!(entries,[:(getproperty(extra,$(QuoteNode(name)))) for name in ExtraNames])
    return :(ComputationDetails($output(($(entries...),))))
end

Base.@constprop :aggressive function _modal_completion(parameters, formulation, selected, diagnostics)
    source_record=parameters.details.data
    source_fields=get(source_record,:formulation_fields,(;))
    modal_descriptions=Engine.completed_formulation(formulation).formulation_fields.all
    modal_descriptions=filter(field -> !isempty(field.meaning),modal_descriptions)
    descriptions=vcat(get(source_fields,:all,NamedTuple[]),modal_descriptions)
    fields=merge(source_fields,
        (all=descriptions,Z=copy(descriptions),Y=copy(descriptions)))
    added=NamedTuple{(:source_gridpoint,:phase_coordinates,:coordinates,:modal,:formulation_fields),
        Tuple{Any,Any,Vector{String},NamedTuple,NamedTuple}}((
         get(parameters.details.data, :gridpoint, nothing),
         get(parameters.details.data, :coordinates, nothing),
         string.(1:size(parameters.Z,1)),
         (identifier=formula_id(selected), requested=NamedTuple(formulation.definition),
            effective=NamedTuple(selected), diagnostics=diagnostics),
         fields))
    return _modal_detail_merge(source_record,added)
end

function compute(::LineCableModelsModal,
        problem::ModalAnalysisProblem{<:LineParameters{T,U,PhaseDomain,Basis}},
        formulation::ModalAnalysisFormulation, execution::ComputationOptions) where {T,U,Basis}
    parameters=problem.parameters
    selected = formulation.formula
    workspace = ModalAnalysisWorkspace(parameters, selected)
    decompose!(allocation_selector(selected), workspace, formula_parameters(selected),
        formulation_options(selected))
    maps = _check_operators(ModalOperators(copy(workspace.Tv), copy(workspace.Ti)), parameters)
    _coordinate_algebra!(workspace, execution)
    roots = _dependent_roots(workspace)
    all(isfinite, roots) || throw(DomainError(roots, "modal roots must be finite"))
    modal = ModalDomain(maps, roots)
    diagnostics = map(copy,workspace.diagnostics)
    retained = _modal_completion(parameters, formulation, selected, diagnostics)
    return LineParameters(modal,
        SeriesImpedance{eltype(workspace.Zm),Basis}(copy(workspace.Zm)),
        ShuntAdmittance{eltype(workspace.Ym),Basis}(copy(workspace.Ym)),
        parameters.f, retained)
end

"""Restore phase coefficients using the retained modal-to-phase bases."""
function transform(::Type{PhaseDomain}, parameters::LineParameters{T,U,D,Basis}) where {T,U,D<:ModalDomain,Basis}
    maps = _check_operators(parameters.domain.operators, parameters)
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
    retained = ComputationDetails(merge(parameters.details.data,
        (coordinates=get(parameters.details.data,:phase_coordinates,nothing),)))
    return LineParameters(PhaseDomain,
        SeriesImpedance{S,Basis}(impedance),
        ShuntAdmittance{S,Basis}(admittance), parameters.f,
        retained)
end

function compute(problem::ModalAnalysisProblem{P}, formulation::ModalAnalysisFormulation;
        options::Union{NamedTuple,ComputationOptions}=ComputationOptions()) where {P}
    return compute(LineCableModelsModal(), problem, formulation; options)
end

function compute(::LineCableModelsModal, problem::ModalAnalysisProblem{P},
        formulation::ModalAnalysisFormulation;
        options::Union{NamedTuple,ComputationOptions}=ComputationOptions()) where {P}
    options = options isa NamedTuple ? ComputationOptions(options) : options
    execution = computation_options(LineCableModelsModal, options)
    result = if execution.data.timing
        measured = Base.@timed compute(LineCableModelsModal(),problem,formulation,execution)
        Engine.retain_gridpoint(measured.value, details(measured.value).data.gridpoint;
            fields=(timing=(wall_seconds=measured.time, bytes=measured.bytes,
                gc_seconds=measured.gctime, compile_seconds=measured.compile_time,
                recompile_seconds=measured.recompile_time),))
    else
        compute(LineCableModelsModal(),problem,formulation,execution)
    end
    execution.data.on_result === nothing || execution.data.on_result(problem, 1, result)
    return result
end

function computation_details(::Type{<:ModalAnalysisFormulation}, result::LineParameters)::ComputationDetails
    return details(result)
end
