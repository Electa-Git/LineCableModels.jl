# Audit sampling coverage independently of the earth formulas. This diagnostic
# records behavior; it is not a production sampling implementation.
# julia --startup-file=no --compiled-modules=existing --project=. test/gauntlet/spectral_sampling_audit.jl
using LineCableModels
using QuadGK
using SpecialFunctions: erfc
using JSON3
using SHA

const E = LineCableModels.Engine
directory = normpath(joinpath(@__DIR__, "..", "..", ".linecablemodels", "qa", "spectral-sampling"))
mkpath(directory)
centre, width, height = 9.0, 1e-4, 1.0
raw(x) = exp(-((x-centre)/width)^2/2)
exact = sqrt(pi/2)*width*exp(-height*centre + (height*width)^2/2)*
    erfc((height*width^2-centre)/(sqrt(2)*width))
results = []
for method in (:quad,:trapz,:cim)
    samples = ComplexF64[]
    kernel(x) = (push!(samples,complex(x)); complex(raw(x)))
    integral = E.SpectralIntegral(Val(:cosine),kernel,
        (height=height,separation=0.0),1.0)
    controls = E.computation_options(E.SpectralIntegral,
        (method=method,options=(rtol=1e-6,atol=1e-12)))
    result = try
        value = E.integrate(controls.method,integral,controls.options,nothing)
        (method=String(method),returned=true,value_real=real(value),value_imag=imag(value),
            absolute_error=abs(value-exact),within_requested_tolerance=abs(value-exact)<=max(1e-12,1e-6*exact))
    catch error
        (method=String(method),returned=false,error=sprint(showerror,error))
    end
    push!(results,merge(result,(evaluations=length(samples),
        nearest_sample_in_widths=minimum(abs.(samples.-centre))/width)))
    println(last(results))
end

# Supplying the feature's location and width resolves it using the same QuadGK
# package. This checks the analytic reference and isolates missing coverage.
breakpoints = [0.0,1.0,centre-12width,centre,centre+12width,16.0,Inf]
resolved,error = quadgk(x->raw(x)*exp(-height*x),breakpoints...;rtol=1e-10,atol=1e-16)
@assert abs(resolved-exact)<1e-14*exact+1e-16
println("Exact narrow-feature integral = ",exact,"; explicit-feature quadrature = ",resolved)

# Scales of the previously observed three-wire, 1 MHz CIM problem.
s = 2pi*1e6*im
epsilon,mu = 8.8541878128e-12,4pi*1e-7
air_branch = sqrt(-s*mu*(s*epsilon))
earth_scale = abs(sqrt(s*mu*(10+s*epsilon)))
cim_samples = 128
scales = (air_branch_real=real(air_branch),earth_scale=earth_scale,
    initial_pencil_step=earth_scale/(cim_samples-1),
    initial_log_amplitude_node=earth_scale/cim_samples,
    initial_log_validation_node=earth_scale/(cim_samples+1),
    note="The initial nonzero nodes in the smallest pencil, amplitude-fit and held-out grids all lie beyond the air transition scale. This is a coverage risk, not proof it is the sole CIM failure cause.")
println(scales)
open(joinpath(directory,"audit.json"),"w") do io
    JSON3.write(io,(source_sha256=bytes2hex(sha256(read(joinpath(@__DIR__,"..","..","src","engine","integration.jl")))),
        narrow_feature=(centre=centre,width=width,height=height,exact=exact,
            explicit_feature_quad=resolved,estimated_error=error,results=results),earth_scales=scales))
end
