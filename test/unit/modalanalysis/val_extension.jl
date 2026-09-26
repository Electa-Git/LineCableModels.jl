@testitem "ModalAnalysis / complete Val formula extension" tags=[:unit] begin
    import LineCableModels.Engine: initialize_buffers, description
    import LineCableModels.Grammar: formulation_options, FormulationOptions
    import LineCableModels.ModalAnalysis: decompose!, Formula
    import LineCableModels: FormulaMethod

    allocations=Ref(0)
    calculations=Ref(0)
    description(::Type{<:Formula{:diagonal_example}};compact=false)=
        compact ? "diagonal example" : "one-mode diagonal example"
    formulation_options(::FormulaMethod{<:Formula{:diagonal_example},typeof(decompose!)}) =
        FormulationOptions()
    function initialize_buffers(::Val{:diagonal_example},::Type{T},input,invariants,
            common) where {T<:Complex}
        allocations[]+=1
        invariants.n==1 || throw(DimensionMismatch("diagonal example requires one mode"))
        return merge(common,(diagonal_product=Vector{T}(undef,invariants.nf),))
    end
    function decompose!(::Val{:diagonal_example},workspace,parameters::NamedTuple,
            options::FormulationOptions)
        calculations[]+=1
        scratch=workspace.buffers.diagonal_product
        for k in eachindex(scratch)
            z=workspace.input.Z[1,1,k]/workspace.input.root_scale
            y=workspace.input.Y[1,1,k]/workspace.input.root_scale
            scratch[k]=z*y
            root=sqrt(scratch[k])
            (real(root)<0 || (iszero(real(root)) && imag(root)<0)) && (root=-root)
            workspace.roots[1,k]=root*workspace.input.root_scale
            workspace.Tv[1,1,k]=one(root)
            workspace.Ti[1,1,k]=one(root)
            workspace.diagnostics.eigen_residual[1,k]=zero(real(root))
            workspace.diagnostics.iterations[1,k]=0
            workspace.diagnostics.converged[1,k]=true
        end
        return workspace
    end

    phase=LineParameters(reshape(ComplexF64[2+im,3+im],1,1,2),
        reshape(ComplexF64[1e-6im,2e-6im],1,1,2),[50.0,100.0];
        details=ComputationDetails((inputs=(system=(line_length=10.0,),),)))
    selected=ModalAnalysisFormulation(Formula(Val(:diagonal_example)))
    modal=compute(ModalAnalysisProblem(phase),selected)
    segment=PropagationParameters(modal)
    @test allocations[]==1 && calculations[]==1
    @test size(gamma(modal))==(1,2)
    @test size(H(segment))==size(gamma(modal))
    @test Tv(modal)==ones(ComplexF64,1,1,2)
    @test Ti(modal)==ones(ComplexF64,1,1,2)
    @test Z(modal)≈Z(phase)
    @test Y(modal)≈Y(phase)
    @test all(details(modal).data.modal.diagnostics.converged)
end
