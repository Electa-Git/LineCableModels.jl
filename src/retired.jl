const LEGACY_FEM_SECTOR_BRANCH = "legacy/fem-sector"
const LEGACY_FEM_SECTOR_COMMIT = "b75dd2723f90a83ec090b20605ea42af57f4a9c3"
const LEGACY_FEM_SECTOR_URL = "https://github.com/Electa-Git/LineCableModels.jl.git"
const LEGACY_FEM_SECTOR_INSTALL = "julia --project=. -e \"using Pkg; " *
                                  "Pkg.add(Pkg.PackageSpec(url=ARGS[1], " *
                                  "rev=ARGS[2]))\" " *
                                  LEGACY_FEM_SECTOR_URL * " " *
                                  LEGACY_FEM_SECTOR_COMMIT
const LEGACY_JSON_COMMIT = "a71bdfe1ac832f27a0c88b1d02596194aac46ec7"

function retired_fem_sector(feature::AbstractString)
    throw(ArgumentError(
        "$feature was removed from LineCableModels. The final FEM and " *
        "sector-support snapshot is branch $LEGACY_FEM_SECTOR_BRANCH at commit " *
        "$LEGACY_FEM_SECTOR_COMMIT.\n\nInstall that snapshot in the current " *
        "Julia project with:\n$LEGACY_FEM_SECTOR_INSTALL",
    ))
end

function retired_legacy_json()
    Base.depwarn(
        "Unversioned LineCableModels JSON is retired. Commit $LEGACY_JSON_COMMIT " *
        "is the last snapshot that can load and migrate this file.",
        :load!,
        force = true
    )
    throw(ArgumentError(
        "legacy LineCableModels JSON is unsupported; use commit " *
        "$LEGACY_JSON_COMMIT to migrate the file",
    ))
end

struct _RetiredFEMNamespace end

function Base.getproperty(::_RetiredFEMNamespace, name::Symbol)
    return retired_fem_sector("Engine.FEM.$name")
end

const _RETIRED_FEM = _RetiredFEMNamespace()

SectorParams(args...; kwargs...) = retired_fem_sector("SectorParams")
Sector(args...; kwargs...) = retired_fem_sector("Sector")
SectorInsulator(args...; kwargs...) = retired_fem_sector("SectorInsulator")
