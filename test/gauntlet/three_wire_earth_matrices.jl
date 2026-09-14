# Manual extension of the earth-only test to three collinear, equally spaced wires.
# julia --startup-file=no --compiled-modules=existing --project=. test/gauntlet/three_wire_earth_matrices.jl
module ThreeWireEarthMatrices

include("earth_matrices_manual.jl")
using .EarthMatricesManual: earth_matrices, registered_xue, GEOMETRY, FREQUENCIES,
    encoded, relative, besseli
using LinearAlgebra
using JSON3
using Printf
using SHA

const OUTPUT = normpath(joinpath(@__DIR__,"..","..",".linecablemodels","qa","three-wire-earth-matrices"))
const MODELS = (:proposal,:xue,:field_average_cf)
const QUANTITIES = (:Ze,:Pe,:Ye)
const ENTRIES = ((1,1),(1,2),(1,3),(2,1),(2,2))
assemble(self,adjacent,distant) = [self adjacent distant; adjacent self adjacent; distant adjacent self]
pack(models) = NamedTuple{MODELS}(Tuple(
    NamedTuple{QUANTITIES}(Tuple(encoded(model[q]) for q in QUANTITIES)) for model in models))

function three_wire_matrices(f,method=:quad;rtol=1e-10,atol=0.0)
    # Reuse the independently audited source kernels at separations 1 m and 2 m.
    # Finalized two-wire P/Y and the full-current two-wire Z are never assembled.
    adjacent = earth_matrices(f,method;rtol,atol)
    distant = earth_matrices(f,method;g=merge(GEOMETRY,(separation=2.0,)),rtol,atol)
    K = assemble(adjacent.K[1,1],adjacent.K[1,2],distant.K[1,2])
    H = assemble(adjacent.H[1,1],adjacent.H[1,2],distant.H[1,2])
    s = complex(0.0,2pi*f)
    shg = inv(GEOMETRY.rho)+s*GEOMETRY.epsilon
    kg = sqrt(s*GEOMETRY.mu*shg)
    A,D = adjacent.A,adjacent.D
    F = 2pi*shg*GEOMETRY.radius*besseli(1,kg*GEOMETRY.radius)/(kg*A)
    identity = Matrix{ComplexF64}(I,3,3)
    L = identity/A-F*K
    function xue_matrix(quantity)
        a,b = adjacent.models.xue[quantity],distant.models.xue[quantity]
        assemble(a[1,1],a[1,2],b[1,2])
    end
    xue_Z,xue_P = xue_matrix(:Ze),xue_matrix(:Pe)
    pair(Z,P) = (Ze=Z,Pe=P,Ye=s*identity/P)
    models = (proposal=pair(K/L,H/L),xue=pair(xue_Z,xue_P),field_average_cf=pair(K/D,H/D))
    proposal = models.proposal
    # Use ordinary transpose for reciprocity checks of complex matrices.
    asymmetry = proposal.Ye-transpose(proposal.Ye)
    commutator_prediction = s*F*(H\(K*H-H*K))/H
    checks = (inverse=maximum(relative(v.Ye*v.Pe,s*identity) for v in models),
        current_Z=relative(proposal.Ze*L,K),
        current_P=relative(proposal.Pe*L,H),
        current_Y=relative(proposal.Ye*H,s*L),
        Z_reciprocity=relative(proposal.Ze,transpose(proposal.Ze)),
        centrosymmetry=maximum(relative(v[q],v[q][3:-1:1,3:-1:1]) for v in models for q in QUANTITIES),
        xue_Y_reciprocity=relative(models.xue.Ye,transpose(models.xue.Ye)),
        Y_commutator=norm(asymmetry-commutator_prediction)/norm(proposal.Ye))
    @assert maximum(checks) < 2e-12
    # This example can expose nonsymmetric path-voltage P and Y. Do not edit them.
    diagnostics = (Y_asymmetry=norm(asymmetry)/norm(proposal.Ye),
        Pe_asymmetry=relative(proposal.Pe,transpose(proposal.Pe)),
        condition_L=cond(L),condition_H=cond(H),condition_Pe=cond(proposal.Pe),
        adjacent_Y_directional=relative(proposal.Ye[1,2],proposal.Ye[2,1]),
        Z_vs_Xue=relative(proposal.Ze,models.xue.Ze),
        Y_vs_Xue=relative(proposal.Ye,models.xue.Ye),
        # A three-wire solve must not simply reproduce a finalized two-wire block.
        Z_vs_two_wire=relative(proposal.Ze[1:2,1:2],adjacent.models.full_current.Ze),
        Y_vs_two_wire=relative(proposal.Ye[1:2,1:2],adjacent.models.full_current.Ye))
    return (;models,K,H,L,A,D,F,checks,diagnostics)
end

function check_xue(f,result)
    adjacent = registered_xue(f)
    distant = registered_xue(f;g=merge(GEOMETRY,(separation=2.0,)))
    make(q) = assemble(adjacent[q][1,1],adjacent[q][1,2],distant[q][1,2])
    Z,P = make(:Ze),make(:Pe)
    Y = complex(0.0,2pi*f)*Matrix{ComplexF64}(I,3,3)/P
    error = maximum(relative(result.models.xue[q],v) for (q,v) in ((:Ze,Z),(:Pe,P),(:Ye,Y)))
    @assert error < 2e-9
    return error
end

function main()
    BLAS.set_num_threads(1)
    mkpath(OUTPUT)
    rows = []
    methods = []
    for f in FREQUENCIES
        result = three_wire_matrices(f)
        xue_error = check_xue(f,result)
        push!(rows,(frequency=f,result,xue_error))
        @printf("%9g Hz: Y12=%+.9e%+.9ej, Y21=%+.9e%+.9ej; directional difference=%.4g%%\n",
            f,real(result.models.proposal.Ye[1,2]),imag(result.models.proposal.Ye[1,2]),
            real(result.models.proposal.Ye[2,1]),imag(result.models.proposal.Ye[2,1]),
            100result.diagnostics.adjacent_Y_directional)
        flush(stdout)
        for method in (:trapz,:cim)
            try
                trial = three_wire_matrices(f,method;rtol=1e-6,atol=method===:cim ? 1e-6 : 0.0)
                errors = [(;model,quantity,
                    matrix_relative=relative(trial.models[model][quantity],result.models[model][quantity]),
                    entry_relative=maximum(relative(trial.models[model][quantity][i,j],
                        result.models[model][quantity][i,j]) for (i,j) in ENTRIES))
                    for model in MODELS for quantity in QUANTITIES]
                push!(methods,(;frequency=f,method,status="ok",errors))
            catch error
                push!(methods,(;frequency=f,method,status="failed",error=sprint(showerror,error)))
            end
        end
    end
    dense = [(frequency=f,models=pack(three_wire_matrices(f).models))
        for f in 10.0 .^ range(-1,6;length=281)]
    data = (geometry=(positions_m=[[0.0,-1.0],[1.0,-1.0],[2.0,-1.0]],
            radius_m=GEOMETRY.radius,earth_resistivity_ohm_m=GEOMETRY.rho,
            relative_permittivity=1,relative_permeability=1),Gamma=0,
        reference="infinite earth depth",units=(Ze="ohm/m",Pe="m/F",Ye="S/m"),
        conventions=(proposal="full manuscript: K/L, H/L, jomega L/H",
            xue="retained underground Xue Ze and Pe; full 3x3 Ye=jomega/Pe",
            field_average_cf="isolated primary-current normalization from the supplied two-wire table"),
        frequencies=FREQUENCIES,source_sha256=bytes2hex(sha256(read(@__FILE__))),
        kernel_source_sha256=bytes2hex(sha256(read(joinpath(@__DIR__,"earth_matrices_manual.jl")))),
        rows=[(frequency=r.frequency,models=pack(r.result.models),
            K=encoded(r.result.K),H=encoded(r.result.H),L=encoded(r.result.L),
            checks=r.result.checks,diagnostics=r.result.diagnostics,registered_xue_error=r.xue_error) for r in rows],
        methods,dense)
    write(joinpath(OUTPUT,"matrices.json"),JSON3.write(data))
    open(joinpath(OUTPUT,"matrices.csv"),"w") do io
        println(io,"frequency_Hz,model,quantity,row,column,real,imag")
        for r in rows, model in MODELS, quantity in QUANTITIES, j in 1:3, i in 1:3
            value = r.result.models[model][quantity][i,j]
            println(io,join((r.frequency,model,quantity,i,j,real(value),imag(value)),','))
        end
    end
    println("Results: ",OUTPUT)
    for entry in methods
        entry.status=="ok" || println(entry)
    end
end

abspath(PROGRAM_FILE)==(@__FILE__) && main()
end
