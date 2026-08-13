using Pkg

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path = ROOT_DIR))
Pkg.instantiate()
