"""Run an explicit PEC boundary reference using a saved bare-wire FEM mesh.

This is an audit-only specialization of the existing quasi-TEM exterior PDE.
It does not use the manuscript's Green kernels and does not modify the backend.
Usage: python3 test/gauntlet/pec_boundary_reference.py PEC_AUDIT_DIRECTORY
"""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np


def encoded(a):
    return dict(shape=list(a.shape), real=a.real.ravel(order="F").tolist(),
                imag=a.imag.ravel(order="F").tolist())


def main(directory):
    directory = directory.resolve()
    finite = json.loads((directory / "fem.json").read_text())
    assert finite["core_material"] == "pec"
    original = Path(finite["details"]["fem"]["run"]["run_directory"])
    settings = json.loads((original / "input/computation.json").read_text())
    assert len(settings["materials"]) == 2
    assert all(m["kind"] == "conductor" for m in settings["materials"])
    root = directory / "pec-boundary-reference"
    root.mkdir(exist_ok=True)
    run = Path(tempfile.mkdtemp(prefix="run-", dir=root))
    source = original / "input/getdp"
    owned = run / "getdp"
    shutil.copytree(source, owned)
    text = (source / "quasi-tem.pro").read_text()
    operator = Path(__file__).parent / "getdp/pec_boundary.pro"
    header = text[:text.index("FunctionSpace {")]
    driver = text[text.index("Macro FEMSetBasisCurrent"):text.index("PostProcessing {")]
    post = '''PostProcessing {
  { Name FEMFields; NameOfFormulation FEM_a_phi_2D; NameOfSystem Sys_FEM;
    PostQuantity {
      { Name ReZ; Value { Term { [Re[Complex[0.,omega[]]*{A}/UnitSource]]; In Terminals; } }}
      { Name ImZ; Value { Term { [Im[Complex[0.,omega[]]*{A}/UnitSource]]; In Terminals; } }}
      { Name ReP; Value { Term { [Re[{Psi}/UnitSource]]; In Terminals; } }}
      { Name ImP; Value { Term { [Im[{Psi}/UnitSource]]; In Terminals; } }}
    }
  }
}
'''
    footer = text[text.index("// Expand output declarations"):]
    (owned / "quasi-tem.pro").write_text(header + operator.read_text() + driver + post + footer)
    materials = (owned / "materials.pro").read_text()
    expected = "DomainC = Region[{ConductorMaterialRegions, Earth, EarthInf}];"
    assert expected in materials
    (owned / "materials.pro").write_text(materials.replace(expected,
        "DomainC = Region[{Earth, EarthInf}];"))
    model = (owned / "model.pro").read_text()
    # The normalized equations take the analytic Gamma -> 0 limit explicitly.
    (owned / "model.pro").write_text(model.replace("GammaQuasiTEMIm = 1.0e-12;",
        "GammaQuasiTEMIm = 0.0;"))
    # inv_gamma is not used by this formulation or its output operation.
    qt = (owned / "quasi-tem.pro").read_text()
    (owned / "quasi-tem.pro").write_text(qt.replace("inv_gamma[] = 1. / gamma_prop[];", ""))
    executable = settings["getdp_provenance"]["path"]
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OPENBLAS_DEFAULT_NUM_THREADS="1",
               OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")

    def solve(plan):
        index = plan["frequency_index"]
        job = run / f"f{index:04d}"
        (job / "raw/jobs").mkdir(parents=True)
        (job / "maps").mkdir()
        (job / "bases.pro").write_text("RequestedBases() = {1,2};\n")
        mesh_name = "model.msh" if index == len(settings["mesh_plans"]) else f"frequency_{index:04d}.msh"
        mesh = original / "mesh" / mesh_name
        assert mesh.exists(), mesh
        command = [executable, str(owned / "model.pro"), "-solve", "LineCableModelsFEMScan",
            "-msh", str(mesh), "-name", str(job / "solver"), "-v", "4",
            "-setstring", "ModelDataPath", str(original / "input/model_data.pro"),
            "-setstring", "RunDirectory", str(job),
            "-setstring", "BasisListPath", str(job / "bases.pro"),
            "-setnumber", "FrequencyIndex", str(index),
            "-setnumber", "FrequencyHz", str(plan["frequency"]),
            "-setnumber", "Val_Rint", str(plan["domain_radius"]),
            "-setnumber", "Val_Rext", str(plan["shell_outer_radius"]),
            "-setnumber", "PlotFieldMaps", "0", "-setnumber", "ReuseFactorization", "1"]
        (job / "command.json").write_text(json.dumps(command, indent=2))
        with (job / "getdp.log").open("w") as log:
            subprocess.run(command, cwd=job, env=env, stdout=log, stderr=subprocess.STDOUT,
                           check=True)
        matrices = []
        for quantity in ("Z", "P"):
            a = np.empty((2, 2), dtype=complex)
            for basis in (1, 2):
                path = job / "raw/jobs" / f"getdp-f{index:04d}-b{basis:04d}-{quantity}.tsv"
                values = np.loadtxt(path, ndmin=2)
                assert values.shape == (2, 6), values.shape
                for fi, f, row, col, re, im in values:
                    assert fi == index and f == plan["frequency"] and col == basis
                    a[int(row)-1, int(col)-1] = complex(re, im)
            assert np.isfinite(a).all()
            matrices.append(a)
        z, p = matrices
        y = np.linalg.solve(p, np.eye(2))
        assert np.all(np.diag(z).real > 0) and np.all(np.diag(z).imag > 0), z
        assert np.linalg.norm(z-z.T)/np.linalg.norm(z) < 1e-8
        print(f"PEC boundary {plan['frequency']:g} Hz: Z11={z[0,0]}, Y11={y[0,0]}", flush=True)
        return z, y, p

    with ThreadPoolExecutor(max_workers=2) as pool:
        values = list(pool.map(solve, settings["mesh_plans"]))
    z, y, p = (np.stack([v[i] for v in values], axis=2) for i in range(3))
    result = dict(frequencies=finite["frequencies"], Z=encoded(z), Y=encoded(y), P=encoded(p),
        core_material="pec", representation="exact PEC boundary; normalized Gamma=0 quasi-TEM",
        geometry_source=str(original), run_directory=str(run),
        operator_sha256=hashlib.sha256(operator.read_bytes()).hexdigest(),
        inversion_residual=max(float(np.linalg.norm(p[:,:,i]@y[:,:,i]-np.eye(2))) for i in range(p.shape[2])))
    (directory / "fem-pec-boundary.json").write_text(json.dumps(result))
    print(directory / "fem-pec-boundary.json")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
