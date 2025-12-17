
using Pkg, Dates

# -------- sources --------
lcm_url = "https://github.com/Electa-Git/LineCableModels.jl"
lcm_rev = "dev/thermal"
getdp_url = "https://github.com/Electa-Git/GetDP.jl"
getdp_rev = "473dece8e9bdb614990f611cac346d10984cfd25"
# -------------------------

Pkg.activate()

# 1) Remove from manifest (works even when not in Project.toml)
for name in ("GetDP", "LineCableModels")
	try
		Pkg.rm(name; mode = Pkg.PKGMODE_MANIFEST)
	catch err
		@warn "Pkg.rm failed (ok if not present)" name err
	end
end

# 2) GC immediately
try
	Pkg.gc(; collect_delay = Day(0))
catch err
	@warn "Pkg.gc failed (non-fatal)" err
end

# 3) Delete depot caches (sources + compiled) for these two packages
jver = "v$(VERSION.major).$(VERSION.minor)"

for depot in DEPOT_PATH
	# Sources / dev checkouts
	for p in (
		joinpath(depot, "packages", "GetDP"),
		joinpath(depot, "packages", "LineCableModels"),
		joinpath(depot, "dev", "GetDP"),
		joinpath(depot, "dev", "LineCableModels"),
	)
		if ispath(p)
			rm(p; force = true, recursive = true)
			println("rm -rf ", p)
		end
	end

	# Compiled caches
	compdir = joinpath(depot, "compiled", jver)
	if isdir(compdir)
		for (root, _, files) in walkdir(compdir)
			for f in files
				if startswith(f, "GetDP") || startswith(f, "LineCableModels")
					fp = joinpath(root, f)
					rm(fp; force = true)
					println("rm ", fp)
				end
			end
		end
	end
end

# 4) Re-add:  pinned GetDP + LineCableModels, then rebuild
Pkg.add(PackageSpec(url = getdp_url, rev = getdp_rev))
Pkg.add(PackageSpec(url = lcm_url, rev = lcm_rev))


Pkg.resolve()
Pkg.instantiate()
Pkg.precompile()

# 5) Verify the resolved GetDP source (truth)
Pkg.status("GetDP")
println(only([i for (u, i) in Pkg.dependencies() if i.name == "GetDP"]).source)
