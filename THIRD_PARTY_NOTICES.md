# Third-party notices

## GetDP

The optional Gmsh/GetDP finite-element backend can download and execute GetDP
3.5.0 from the package's lazy `getdp` artifact. GetDP is copyright (C)
1997–2022 P. Dular and C. Geuzaine, University of Liege, and is distributed
under the GNU General Public License, version 2 or later.

- Project: https://getdp.info/
- Source: https://getdp.info/src/getdp-3.5.0-source.tgz
- Source SHA-256: `d6814dc3f81431f1db30b3d5318553efab616d7ea53b352a2c2d0640d130a328`
- License: https://getdp.info/doc/texinfo/getdp.html#License

The upstream binary archives bound by `Artifacts.toml` include `LICENSE.txt`,
`CREDITS.txt`, and `README.txt`. LineCableModels invokes GetDP as an external
program; GetDP is not incorporated into the LineCableModels library.

Maintainers can revalidate every supported archive, tree hash, executable and
license file with `julia dev/verify_getdp_artifact.jl`.
