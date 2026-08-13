# Third-party notices

## Gmsh

The optional FEM extension uses the registered Gmsh.jl package. Its Julia
wrapper is MIT-licensed; the Gmsh library and executable are distributed under
the GNU General Public License, version 2 or later, with Gmsh's linking
exception. Gmsh is resolved by Julia's package manager and is not vendored in
this repository. See <https://gmsh.info/> for source and licensing.

## GetDP frontend snapshot

The private Julia frontend under `ext/fem/getdp_frontend/` is derived from
Electa-Git/GetDP.jl commit
`b1d91b0d8974ea642b772462edcf6e26299fdf0a`. It is distributed under the BSD
3-Clause License included in that directory. `NOTICE.md` records the upstream
tree and the mechanical compatibility changes.

## GetDP executable

The optional FEM integration can invoke a separately installed GetDP
executable. GetDP 3.5.0 is distributed by its authors under the GNU General
Public License, version 2 or later. LineCableModels does not include, download,
or redistribute GetDP binaries or source code. Users provide the executable
through `GETDP_EXECUTABLE` or `PATH` and remain responsible for its license
terms. See <https://www.getdp.info/> for GetDP downloads, source, and licensing.
