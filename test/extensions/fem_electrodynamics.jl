@testitem "Gmsh FEM / basic buried wires / native result transport and extraction" tags=[:extension,:integration,:fem_numerical] setup=[GauntletSupport] begin
    using .GauntletSupport.Gauntlet
    using Gmsh, LinearAlgebra
    reductions=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false)
    selected=Formulation(:LineCableModelsFEM;options=(;reductions...,physics=:quasi_tem))
    controls=(mesh_policy=:remesh,keep_run_directory=true,frequency_workers=1,
        solver_threads=1,gmsh_verbosity=0,getdp_verbosity=0,output_basis=:pul,trace=true)
    # The existing cases exercise actual native solves. Returned quantities must
    # preserve terminal/frequency identity and the declared extraction algebra.
    # Physical agreement with other formulations belongs to Gauntlet studies.
    for case_id in (:two_bare_wires,:two_insulated_wires)
        @testset "$case_id" begin
            model=load_case(case_id;variation=ExactOverrides(frequencies=[50.0,1e4]))
            problem=model.problem
            elapsed=@elapsed actual=compute(problem,selected;options=controls)
            record=details(actual).fem
            @test frequencies(actual)==problem.frequencies
            @test details(actual).coordinates==model.port_order
            @test basis(actual)===:pul
            @test domain(actual)===PhaseDomain
            @test size(Z(actual))==size(Y(actual))==(2,2,2)
            @test all(isfinite,Z(actual)) && all(isfinite,Y(actual))
            @test record.inputs.options.physics===Symbol("quasi-tem")
            @test record.inputs.execution.domain_skin_depths==2.0
            @test record.primitive.phase_map==record.reduced_phase_map==problem.system.connection_order
            @test record.run.completed_frequencies==2
            @test record.run.completed_columns==4
            @test Z(actual)==record.primitive.Z_primitive
            for k in eachindex(problem.frequencies)
                @test Y(actual)[:,:,k] ≈ inv(record.primitive.P_primitive[:,:,k]) rtol=1e-12
            end
            # Bind every returned primitive entry to the native output protocol.
            for (file,expected) in (("Z.tsv",record.primitive.Z_primitive),
                    ("P.tsv",record.primitive.P_primitive))
                lines=readlines(joinpath(record.run.run_directory,"raw",file))
                @test length(lines)==1+length(expected)
                seen=Set{NTuple{3,Int}}()
                for line in lines[2:end]
                    columns=split(line,'\t')
                    @test length(columns)==6
                    k,i,j=parse.(Int,columns[[1,3,4]])
                    @test parse(Float64,columns[2])==problem.frequencies[k]
                    @test (i,j,k) ∉ seen
                    push!(seen,(i,j,k))
                    @test complex(parse(Float64,columns[5]),parse(Float64,columns[6]))==expected[i,j,k]
                end
                @test length(seen)==length(expected)
            end
            println(case_id,": FEM ",elapsed," s; retained run ",record.run.run_directory)
        end
    end
end
