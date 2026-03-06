using Revise
include("CableBuilder.jl")
using .CableBuilder
using BenchmarkTools

ms_cu = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393, 0.0, 0.0, 0.0, 0.0)
ms_vac = Material(Inf, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

parts = (
	Conductor.SolidCore(:core, ms_cu; r = 0.02315),
	Conductor.Tubular(:sheath, ms_cu; t = 0.005),
	Conductor.Pipe(:pipe, ms_cu; t = 0.005, filler = ms_vac, offset = 0.01),
	Conductor.Wires(:sheath, ms_cu; r_w = 0.005, n_w = 3, lay_r = 0.01),
)

# Flawless compilation
cds = CableDesignSpec(parts)

ast = first(cds)
println(ast.payload)


# # 1. Define the Blueprint (The User DSL)
# # Notice we mix deterministic scalars, a relative sweep, and an absolute sweep
# copper_spec = MaterialSpec(
# 	rho = Grid(1.68e-8:0.01e-8:1.72e-8, 2.0), # 2% manufacturing tolerance
# 	eps_r = 1.0,                              # Auto-promoted to DeterministicGrid
# 	mu_r = 1.0,
# 	T0 = Grid(20.0:10.0:90.0, AbsoluteError(2.0)), # ± 2.0°C absolute sensor error
# )

# @show typeof(copper_spec)
# @code_warntype rand(copper_spec)

# @show @allocated rand(copper_spec)
# @btime @allocated rand($copper_spec)


# # # ---------------------------------------------------------
# # # Benchmark 1: The Array Comprehension (User Snippet)
# # # ---------------------------------------------------------
# # println("--- Benchmarking 1000 Monte Carlo Realizations (Array Allocation) ---")
# # # We interpolate $copper_spec to avoid global scope penalty in the benchmark
# # @btime mc_materials = [rand($copper_spec) for _ in 1:1000]

# # # ---------------------------------------------------------
# # # Benchmark 2: The Bare-Metal Loop (Zero-Allocation Target)
# # # ---------------------------------------------------------
# # # This tests the pure speed of your ntuple/Grid/rand architecture 
# # # without the overhead of Julia allocating a 1000-element Vector
# # function pure_monte_carlo_loop(spec, N)
# # 	# We just draw the sample and discard it to test the engine's raw speed
# # 	for _ in 1:N
# # 		m = rand(spec)
# # 	end
# # 	return nothing
# # end

# # println("\n--- Benchmarking 1000 Monte Carlo Realizations (Engine Only) ---")
# # @btime pure_monte_carlo_loop($copper_spec, 1000)


# # ---------------------------------------------------------
# # Test 1: The Deterministic Single Build
# # ---------------------------------------------------------
# println("--- Benchmarking Single Deterministic Build ---")
# # We use $ to interpolate the variable into the macro so it doesn't benchmark global scope lookup
# @btime first($cds)

# # ---------------------------------------------------------
# # Test 2: The Combinatorial Iterator (The Real Test)
# # ---------------------------------------------------------
# grid_parts = (
# 	Conductor.Solid(:core, ms_cu; r = Grid([0.02, 0.025, 0.03])), # 3 variations
# 	Conductor.Tubular(:sheath, ms_cu; t = Grid([0.004, 0.005, 0.006, 0.007])), # 4 variations
# 	Conductor.Pipe(
# 		:pipe,
# 		ms_cu;
# 		t = Grid([0.004, 0.005, 0.006, 0.007]),
# 		filler = ms_vac,
# 		offset = Grid([0.004, 0.005, 0.006, 0.007]),
# 	), # 4 variations
# )
# grid_cds = CableDesignSpec(grid_parts)
# println("\n--- Benchmarking 1-Design Materialization from Grids ---")
# @btime first($grid_cds)

# # A function barrier to test the loop exactly how your solver will use it
# function exhaust_generator(spec)
# 	count = 0
# 	for design in spec
# 		# The design is materialized here. 
# 		# If the compiler is happy, this loop will allocate ZERO memory.
# 		count += 1
# 	end
# 	return count
# end

# println("\n--- Benchmarking Combinatorial Sweep ---")
# @btime exhaust_generator($grid_cds)

# function alloc_per_design(spec, n)
# 	s = Iterators.take(spec, n)
# 	a = @allocated for x in s
# 		nothing
# 	end
# 	return a / n
# end
# println("\n--- Alloc per design ---")
# @show alloc_per_design(grid_cds, length(grid_cds))
