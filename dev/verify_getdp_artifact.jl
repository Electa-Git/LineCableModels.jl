using Downloads
using Pkg
using SHA
using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const ARTIFACTS = TOML.parsefile(joinpath(ROOT, "Artifacts.toml"))["getdp"]
const REQUIRED_FILES = Dict(
    "linux" => (
        "getdp-3.5.0-Linux64/bin/getdp",
        "getdp-3.5.0-Linux64/share/doc/getdp/LICENSE.txt",
    ),
    "macos" => (
        "getdp-3.5.0-MacOSX/bin/getdp",
        "getdp-3.5.0-MacOSX/share/doc/getdp/LICENSE.txt",
    ),
)

for binding in ARTIFACTS
    download = only(binding["download"])
    mktempdir() do directory
        archive = Downloads.download(
            download["url"], joinpath(directory, basename(download["url"])),
        )
        digest = bytes2hex(open(sha256, archive))
        digest == download["sha256"] || error(
            "archive checksum mismatch for $(binding["os"]): $digest",
        )
        unpacked = joinpath(directory, "unpacked")
        Pkg.PlatformEngines.unpack(archive, unpacked)
        tree = bytes2hex(Pkg.GitTools.tree_hash(unpacked))
        tree == binding["git-tree-sha1"] || error(
            "artifact tree mismatch for $(binding["os"]): $tree",
        )
        all(isfile(joinpath(unpacked, path)) for path in REQUIRED_FILES[binding["os"]]) ||
            error("GetDP executable or license is missing for $(binding["os"])")
        println(binding["os"], " ", binding["arch"], ": verified ", tree)
    end
end
