const LEGACY_FEM_SECTOR_BRANCH = "legacy/fem-sector"
const LEGACY_FEM_SECTOR_COMMIT = "b75dd2723f90a83ec090b20605ea42af57f4a9c3"
const LEGACY_FEM_SECTOR_URL = "https://github.com/Electa-Git/LineCableModels.jl.git"
const LEGACY_FEM_SECTOR_INSTALL = "julia --project=. -e \"using Pkg; " *
                                  "Pkg.add(Pkg.PackageSpec(url=ARGS[1], " *
                                  "rev=ARGS[2]))\" " *
                                  LEGACY_FEM_SECTOR_URL * " " *
                                  LEGACY_FEM_SECTOR_COMMIT

function retired_fem_sector(feature::AbstractString)
    throw(
        ArgumentError(
        "$feature was removed from LineCableModels. " *
        "The final FEM and sector-support snapshot is branch " *
        "$(LEGACY_FEM_SECTOR_BRANCH) at commit " *
        "$(LEGACY_FEM_SECTOR_COMMIT).\n\n" *
        "Install that snapshot in the current Julia project with:\n" *
        LEGACY_FEM_SECTOR_INSTALL,
    ),
    )
end
