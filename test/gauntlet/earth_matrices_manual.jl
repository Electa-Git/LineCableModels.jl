# Manual earth-matrix comparison; no cable assembly or FEM dependencies.
# julia --startup-file=no --compiled-modules=existing --project=. test/gauntlet/earth_matrices_manual.jl
module EarthMatricesManual

using LineCableModels
using LinearAlgebra
using SpecialFunctions: besseli, besselk
using JSON3
using Printf
using SHA
const E = LineCableModels.Engine
const FREQUENCIES = 10.0 .^ (-1:6)
const GEOMETRY = (radius=0.0425, depth=1.0, separation=1.0, rho=0.1,
    epsilon=8.8541878128e-12, mu=4pi*1e-7)
const OUTPUT = joinpath(@__DIR__, "..", "..", ".linecablemodels", "qa", "earth-matrices-manual")
const MODELS = (:field_average, :point_field, :xue, :full_current)

relative(a,b) = norm(a-b)/max(norm(b),floatmin(Float64))
encoded(a) = (shape=collect(size(a)), real=vec(real.(a)), imag=vec(imag.(a)))
equal_pair(self,mutual) = ComplexF64[self mutual; mutual self]

function spectral(f, y, kind, method, g; rtol=1e-10, atol=0.0)
    s = complex(0.0,2pi*f)
    sh0, shg = s*g.epsilon, inv(g.rho)+s*g.epsilon
    k02, kg2 = s*g.mu*sh0, s*g.mu*shg
    kernel = lambda -> begin
        a0, ag = sqrt(lambda^2+k02), sqrt(lambda^2+kg2)
        # The SpectralIntegral weight supplies exp(-2h*lambda) cos(y*lambda).
        decay = exp(-2g.depth*kg2/(ag+lambda))
        if kind === :Q
            decay/(a0+ag)
        elseif kind === :MU
            decay*shg*a0/(ag*(shg*a0+sh0*ag))
        else
            error("unknown kernel $kind")
        end
    end
    integral = E.SpectralIntegral(Val(:cosine),kernel,
        (height=2g.depth,separation=y),abs(sqrt(kg2));
        angle=min(pi/4,atan(g.depth/max(y,eps()))))
    controls = E.computation_options(E.SpectralIntegral,
        (method=method, options=(rtol=rtol,atol=atol)))
    return E.integrate(controls.method,integral,controls.options,nothing)
end

function earth_matrices(f, method=:quad; g=GEOMETRY, rtol=1e-10, atol=0.0)
    s = complex(0.0,2pi*f)
    shg = inv(g.rho)+s*g.epsilon
    kg = sqrt(s*g.mu*shg)
    x = kg*g.radius
    A = besseli(0,x)
    D = x*besselk(1,x) # Isolated primary-current map; Cf = 1/D.
    F = 2pi*shg*g.radius*besseli(1,x)/(kg*A)
    integral(y,kind) = spectral(f,y,kind,method,g;rtol,atol)
    # Centre sampling for regular self fields versus Xue's horizontal radius sampling.
    Q0, Qr, Qm = (integral(y,:Q) for y in (0.0,g.radius,g.separation))
    U0, Ur, Um = (integral(y,:MU) for y in (0.0,g.radius,g.separation))
    direct_self = besselk(0,x)
    image_self = besselk(0,kg*2g.depth)
    shifted_image_self = besselk(0,kg*hypot(g.radius,2g.depth))
    mutual_difference = besselk(0,kg*g.separation)-
                        besselk(0,kg*hypot(g.separation,2g.depth))
    K = s*g.mu/(2pi)*equal_pair(
        direct_self+A*(-image_self+2Q0), A*(mutual_difference+2Qm))
    H = s/(2pi*shg)*equal_pair(
        direct_self+A*(-image_self+2U0), A*(mutual_difference+2Um))
    Kpoint = s*g.mu/(2pi)*equal_pair(
        direct_self-shifted_image_self+2Qr, mutual_difference+2Qm)
    Hpoint = s/(2pi*shg)*equal_pair(
        direct_self-shifted_image_self+2Ur, mutual_difference+2Um)
    identity = Matrix{ComplexF64}(I,2,2)
    L = identity/A-F*K # Full manuscript current map at Gamma=0.
    # Pe has units m/F, so Ye = j*omega / Pe, using a full matrix solve.
    pair(Z,P) = (Ze=Z, Pe=P, Ye=s*identity/P)
    models = (field_average=pair(K/D,H/D), point_field=pair(Kpoint/D,Hpoint/D),
        xue=pair(Kpoint,Hpoint), full_current=pair(K/L,H/L))
    checks = (inverse=maximum(relative(v.Ye*v.Pe,s*identity) for v in models),
        point_Z_scaling=relative(models.point_field.Ze,models.xue.Ze/D),
        point_P_scaling=relative(models.point_field.Pe,models.xue.Pe/D),
        point_Y_scaling=relative(models.point_field.Ye,D*models.xue.Ye),
        mutual_average_Z=relative(models.field_average.Ze[1,2],A/D*models.xue.Ze[1,2]),
        mutual_average_P=relative(models.field_average.Pe[1,2],A/D*models.xue.Pe[1,2]),
        current_Z=relative(models.full_current.Ze*L,K),
        current_P=relative(models.full_current.Pe*L,H))
    @assert maximum(checks) < 5e-13
    return (;models,checks,K,H,L,A,D,x)
end

# Independently evaluate the registered, unsimplified S11/S12/S13 equations.
# The former underground defaults are retained explicitly as :Xue2018.
function registered_xue(f; g=GEOMETRY)
    s = complex(0.0,2pi*f)
    matrices = map((E.EarthImpedance,E.EarthAdmittance)) do owner
        formula = owner.Formula(:Xue2018;
            options=(integration=(method=:quad,options=(rtol=1e-11,)),))
        entries = map(1:2) do column
            self = column == 1
            pair = E.EarthPair(1,column,(-g.depth,-g.depth),
                self ? 0.0 : g.separation,(2,2);radius=self ? g.radius : nothing)
            formula([Inf,g.rho],fill(g.epsilon,2),fill(g.mu,2),s,pair)()
        end
        equal_pair(entries...)
    end
    Ze,Pe = matrices
    return (;Ze,Pe,Ye=s*Matrix{ComplexF64}(I,2,2)/Pe)
end

function compare_expected(rows)
    comparisons = []
    lines = readlines(joinpath(@__DIR__,"earth_matrices_expected.tsv"))
    for line in lines[2:end]
        fields = split(line,'\t')
        frequency = parse(Float64,fields[1])
        model,quantity = Symbol(fields[2]),Symbol(fields[3])
        column = parse(Int,fields[4])
        expected = complex(parse(Float64,fields[5]),parse(Float64,fields[6]))
        record = only(filter(row -> row.frequency==frequency,rows))
        actual = getproperty(getproperty(record.result.models,model),quantity)[1,column]
        push!(comparisons,(;frequency,model,quantity,column,
            expected_real=real(expected),expected_imag=imag(expected),
            computed_real=real(actual),computed_imag=imag(actual),
            relative_error=relative(actual,expected),
            real_relative_error=relative(real(actual),real(expected)),
            imag_relative_error=relative(imag(actual),imag(expected))))
    end
    return comparisons
end

function main()
    BLAS.set_num_threads(1)
    mkpath(OUTPUT)
    rows = []
    methods = []
    for f in FREQUENCIES
        result = earth_matrices(f)
        registered = registered_xue(f)
        check = maximum(relative(getproperty(result.models.xue,k),getproperty(registered,k))
            for k in (:Ze,:Pe,:Ye))
        @assert check < 2e-9
        push!(rows,(frequency=f,result,registered_xue_error=check))
        @printf("%9g Hz: registered Xue agreement %.3e; YP residual %.3e\n",f,check,result.checks.inverse)
        flush(stdout)
        for method in (:trapz,:cim)
            try
                trial = earth_matrices(f,method;rtol=1e-6,atol=method===:cim ? 1e-6 : 0.0)
                errors = [(;model,quantity,
                    matrix_relative=relative(getproperty(trial.models[model],quantity),
                        getproperty(result.models[model],quantity)),
                    mutual_relative=relative(getproperty(trial.models[model],quantity)[1,2],
                        getproperty(result.models[model],quantity)[1,2]))
                    for model in MODELS for quantity in (:Ze,:Pe,:Ye)]
                push!(methods,(;frequency=f,method,status="ok",errors))
            catch exception
                push!(methods,(;frequency=f,method,status="failed",error=sprint(showerror,exception)))
            end
        end
    end
    comparisons = compare_expected(rows)
    @printf("Supplied table: %d complex entries; maximum relative difference %.3e\n",
        length(comparisons),maximum(row.relative_error for row in comparisons))
    @assert maximum(row.relative_error for row in comparisons) < 5e-8
    @assert maximum(max(row.real_relative_error,row.imag_relative_error) for row in comparisons) < 2e-6
    dense_frequencies = 10.0 .^ range(-1,6;length=281)
    serialize_models(models) = NamedTuple{MODELS}(Tuple(
        (Ze=encoded(m.Ze),Pe=encoded(m.Pe),Ye=encoded(m.Ye)) for m in models))
    dense = [(frequency=f,models=serialize_models(earth_matrices(f).models)) for f in dense_frequencies]
    data = (geometry=GEOMETRY,Gamma=0,reference="infinite earth depth",units=(Ze="ohm/m",Pe="m/F",Ye="S/m"),
        source_sha256=bytes2hex(sha256(read(@__FILE__))),
        convention=(field_average="circumferential trace / isolated primary current (Cf)",
            point_field="Xue point trace / isolated primary current (Cf)",
            xue="registered underground Xue, algebraically reduced",
            full_current="supplied manuscript: complete L current map"),
        frequencies=FREQUENCIES,
        rows=[(frequency=row.frequency,models=serialize_models(row.result.models),
            K=encoded(row.result.K),H=encoded(row.result.H),L=encoded(row.result.L),
            A=encoded([row.result.A]),D=encoded([row.result.D]),x=encoded([row.result.x]),
            checks=row.result.checks,registered_xue_error=row.registered_xue_error) for row in rows],
        methods,comparisons,dense)
    write(joinpath(OUTPUT,"matrices.json"),JSON3.write(data))
    open(joinpath(OUTPUT,"matrices.csv"),"w") do io
        println(io,"frequency_Hz,model,quantity,row,column,real,imag")
        for row in rows, model in MODELS, quantity in (:Ze,:Pe,:Ye), j in 1:2, i in 1:2
            value = getproperty(row.result.models[model],quantity)[i,j]
            println(io,join((row.frequency,model,quantity,i,j,real(value),imag(value)),','))
        end
    end
    open(joinpath(OUTPUT,"expected-comparison.csv"),"w") do io
        println(io,join(keys(first(comparisons)),','))
        for row in comparisons
            println(io,join(values(row),','))
        end
    end
    println("Results: ",OUTPUT)
    for method in methods
        method.status=="ok" || println(method)
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
end
