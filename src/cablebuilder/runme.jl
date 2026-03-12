using Revise
include("CableBuilder.jl")
using .CableBuilder
import .CableBuilder: CableDesignSpec, MaterialSpec

using BenchmarkTools

# # 1. Define a dummy stochastic Material
# # 2 variations of resistivity
# mat = MaterialSpec(
# 	rho = Grid([1.6e-8, 1.7e-8]),
# 	eps_r = 1.0, mu_r = 1.0, T0 = 20.0, alpha = 0.0,
# 	rho_thermal = 0.0, theta_max = 90.0, tan_delta = 0.0, sigma_solar = 0.0,
# )

# # 2. Define the Solid Core wrapping the Circle primitive
# # 5 variations of radius
# core = Conductor.Solid(:main_core, mat; r = Grid([0.01, 0.02, 0.03, 0.04, 0.05]))

# # 3. Wrap it in the Design Blueprint
# spec = CableDesignSpec((core,))

# # 4. The strict function barrier
# # This forces the compiler to optimize the loop exactly as it would in production
# function exhaust_generator(s)
# 	count = 0
# 	for design in s
# 		# 'design' is fully materialized right here!
# 		# If the primitives, builders, or tuple recursion leak memory, 
# 		# it will pile up on the heap inside this loop.
# 		count += 1
# 	end
# 	return count
# end

# # Warmup to compile
# exhaust_generator(spec)

# # The moment of truth
# println("--- Combinatorial Allocation Test ---")
# @btime exhaust_generator($spec)

ms_cu = Material(1.7241e-8, 1.0, 0.999994, 20.0, 0.00393, 0.0, 0.0, 0.0, 0.0)
ms_vac = Material(Inf, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

parts = (
	Conductor.Solid(:core, ms_cu; r = 0.02315),
	Conductor.Tubular(:sheath, ms_cu; t = 0.005),
	Conductor.Pipe(:pipe, ms_cu; t = 0.005, filler = ms_vac, offset = 0.01),
	# Conductor.Stranded(:sheath, ms_cu; r_w = 0.005, n_w = 3, lay_r = 0.01),
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


# # ---------------------------------------------------------
# # Benchmark 1: The Array Comprehension (User Snippet)
# # ---------------------------------------------------------
# println("--- Benchmarking 100000 Monte Carlo Realizations (Array Allocation) ---")
# @btime mc_materials = [rand($copper_spec) for _ in 1:100000]

# # @allocated(rand(copper_spec)) = 80
# #  27.698 ns (0 allocations: 0 bytes)
# # --- Benchmarking 100000 Monte Carlo Realizations (Array Allocation) ---
# #  2.084 ms (3 allocations: 6.87 MiB)
# #

# # ---------------------------------------------------------
# # Benchmark 2: The Bare-Metal Loop (Zero-Allocation Target)
# # ---------------------------------------------------------
# # This tests the pure speed of your ntuple/Grid/rand architecture 
# # without the overhead of Julia allocating a 1000-element Vector
# function pure_monte_carlo_loop(spec, N)
# 	# We just draw the sample and discard it to test the engine's raw speed
# 	for _ in 1:N
# 		m = rand(spec)
# 	end
# 	return nothing
# end

# println("\n--- Benchmarking 1000 Monte Carlo Realizations (Engine Only) ---")
# @btime pure_monte_carlo_loop($copper_spec, 1000)


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
