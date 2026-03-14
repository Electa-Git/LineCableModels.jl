# diagnose_leaks.jl
using BenchmarkTools
using InteractiveUtils


include("CableBuilder.jl")
using .CableBuilder
import .CableBuilder: CableDesignSpec

println("======================================================")
println("              SROA DEATH INVESTIGATION                ")
println("======================================================\n")

println("--- 1. Setup Dummy Blueprint ---")
mat = Material(
	rho = Grid((1.6e-8, 1.7e-8)),
	eps_r = 1.0, mu_r = 1.0, T0 = 20.0, alpha = 0.0,
	rho_thermal = 0.0, theta_max = 90.0, tan_delta = 0.0, sigma_solar = 0.0,
)

# 5 variations of radius * 2 variations of rho = 10 iterations
core = Conductor.Solid(:core, mat; r = Grid((0.01, 0.02, 0.03, 0.04, 0.05)))
spec = CableDesignSpec(core)

println("\n--- 2. Surgical Allocation Probes ---")
# ROAST: We break the loop apart because `@btime` hides exactly which step is shitting the bed.
function exhaust_generator_diag(s)
	# 1. Setup phase
	alloc_iter = @allocated iter = Base.Iterators.ProductIterator(s.grids)
	println("[Probe] Iterator setup allocs: ", alloc_iter, " bytes")

	iter = Base.Iterators.ProductIterator(s.grids)
	next_val = iterate(iter)

	count = 0
	while next_val !== nothing
		args, state = next_val

		# 2. Construction phase
		# ROAST: If this prints > 0, your abstract Tuple payloads are forcing dynamic dispatch.
		a_build = @allocated design = CableDesign(args...)

		# 3. Advancement phase
		# ROAST: If this prints > 0, your Grid measurement structs are type-unstable.
		a_adv = @allocated next_val = iterate(iter, state)

		if count == 0
			println("[Probe] Loop 1 - CableDesign(args...) allocs: ", a_build, " bytes")
			println("[Probe] Loop 1 - iterate(iter, state) allocs: ", a_adv, " bytes")
		end

		count += 1
	end
	return count
end

exhaust_generator_diag(spec)

println("\n--- 3. Type Stability X-Ray (@code_warntype) ---")
iter = Base.Iterators.ProductIterator(spec.grids)
next_val = iterate(iter)
args, state = next_val

println("[X-Ray] Type of payload args tuple:")
println(typeof(args))
println(
	"\n[X-Ray] If the type above contains `Any` or lacks strict parameters, you already lost.",
)

println("\n------------------------------------------------------")
println("X-Ray: CableDesign(args...)")
println("------------------------------------------------------")
# ROAST: Look for red `Any` or yellow `Union` in the output below. That is the exact variable the compiler boxed.
@code_warntype CableDesign(args...)

println("\n------------------------------------------------------")
println("X-Ray: iterate(iter, state)")
println("------------------------------------------------------")
@code_warntype iterate(iter, state)

println("\n--- 4. Baseline Benchmark ---")
function exhaust_generator(s)
	count = 0
	for design in s
		count += 1
	end
	return count
end

# Warmup
exhaust_generator(spec)

println("Running @btime exhaust_generator(spec)...")
@btime exhaust_generator($spec)