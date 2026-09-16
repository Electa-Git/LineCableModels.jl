"""Current launcher and install contracts, using disposable applications.

Python is the existing test mechanism. These checks establish process and file
ownership, not a language decision or a historical command snapshot.
"""
import json
import os
from pathlib import Path
import selectors
import shutil
import signal
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]


class LauncherContracts(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory(prefix="lcm-contract-")
        self.addCleanup(directory.cleanup)
        self.base = Path(directory.name)
        self.checkout = self.base / "candidate checkout"
        self.checkout.mkdir()
        shutil.copytree(ROOT / "cli", self.checkout / "cli")
        self.version = "0.2.0-validation"
        (self.checkout / "Project.toml").write_text(f'version = "{self.version}"\n')
        self.cwd = self.base / "caller directory"
        self.cwd.mkdir()
        self.bin = self.base / "bin directory"
        self.apps = self.base / "application links"
        self.record = self.base / "invocation.json"
        self.env = dict(os.environ, LCM_INSTALL_DIR=str(self.bin),
                        LCM_APPLICATIONS_DIR=str(self.apps), LCM_RECORD=str(self.record))
        self.launcher = self.checkout / "cli/lcm"
        self.installer = self.checkout / "cli/install.sh"

    def run_command(self, *args, executable=None):
        return subprocess.run([str(executable or self.launcher), *args],
                              cwd=self.cwd, env=self.env, text=True,
                              capture_output=True, timeout=10)

    def app(self, name, *, checkout=None, status=0):
        location = (checkout or self.checkout) / name
        location.mkdir(parents=True)
        script = location / "lcm"
        script.write_text("#!/usr/bin/env python3\n"
                          "import json,os,sys\n"
                          "from pathlib import Path\n"
                          "Path(os.environ['LCM_RECORD']).write_text(json.dumps({"
                          "'cwd':os.getcwd(),'argv':sys.argv,'pid':os.getpid()}))\n"
                          f"sys.exit({status})\n")
        script.chmod(0o755)
        return location

    def test_argument_boundaries_caller_directory_and_exit(self):
        target = self.app("study", status=23)
        self.bin.mkdir()
        (self.bin / "inner").symlink_to(self.launcher)
        public = self.bin / "lcm"
        public.symlink_to("inner")
        args = ["study", "run", "case with spaces", "", "line\nbreak", "$(touch side-effect)", "a'b\"c"]
        result = self.run_command(*args, executable=public)
        self.assertEqual(result.returncode, 23, result.stderr)
        record = json.loads(self.record.read_text())
        self.assertEqual(record["cwd"], str(self.cwd))
        self.assertEqual(record["argv"], [str(target / "lcm"), *args])
        self.assertFalse((self.cwd / "side-effect").exists())

    def test_inspection_and_rejected_commands_do_not_launch(self):
        self.app("study")
        for args in [(), ("--help",), ("-h",), ("help",), ("--paths",)]:
            result = self.run_command(*args)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("study", result.stdout)
            self.assertFalse(self.record.exists())
        self.assertEqual(self.run_command("--version").stdout, f"lcm {self.version}\n")
        for args in [("absent",), ("../study",), ("",), ("cli",), ("--unknown",),
                     ("help", "study", "extra"), ("--version", "extra")]:
            result = self.run_command(*args)
            self.assertEqual(result.returncode, 2)
            self.assertTrue(result.stderr)
            self.assertFalse(self.record.exists())

    def test_application_ownership_and_namespace_forwarding(self):
        outside = self.app("playground", checkout=self.base / "external")
        result = self.run_command("--application", "playground", str(outside), executable=self.installer)
        self.assertEqual(result.returncode, 0, result.stderr)
        for command in ("playground", "runtime", "worker", "presentation", "nats", "container", "demo"):
            self.assertEqual(self.run_command(command, "inspect").returncode, 0)
            self.assertEqual(json.loads(self.record.read_text())["argv"],
                             [str(self.apps / "playground/lcm"), command, "inspect"])
        local = self.app("playground")
        self.assertEqual(self.run_command("help", "playground").returncode, 0)
        self.assertEqual(json.loads(self.record.read_text())["argv"],
                         [str(local / "lcm"), "playground", "--help"])
        paths = self.run_command("--paths").stdout
        self.assertIn(str(local), paths)
        self.assertNotIn(str(outside), paths)

    def test_gauntlet_selects_its_project_and_consumes_namespace(self):
        target = self.checkout / "gauntlet"
        target.mkdir()
        shutil.copy2(ROOT / "gauntlet/lcm", target / "lcm")
        julia = self.app("record-julia", checkout=self.base) / "lcm"
        self.env["LCM_JULIA"] = str(julia)
        self.assertEqual(self.run_command("gauntlet", "inspect", "local case.jl").returncode, 0)
        self.assertEqual(json.loads(self.record.read_text())["argv"], [str(julia),
            f"--project={target}", "--startup-file=no", str(target / "cli.jl"), "inspect", "local case.jl"])

    def test_install_preflight_and_file_preservation(self):
        app = self.app("study")
        self.bin.mkdir()
        public = self.bin / "lcm"
        public.write_bytes(b"user-owned\x00content")
        result = self.run_command("--application", "study", str(app), executable=self.installer)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(public.read_bytes(), b"user-owned\x00content")
        self.assertFalse(self.apps.exists())
        public.unlink()
        public.symlink_to(app / "lcm")
        invalid = [("--bin-dir",), ("--application", "study"),
            ("--application", "runtime", str(app)), ("--application", "../study", str(app)),
            ("--application", "study", str(app), "--application", "study", str(app)),
            ("--application", "study", str(app), "--application", "missing", str(self.base))]
        for args in invalid:
            self.assertNotEqual(self.run_command(*args, executable=self.installer).returncode, 0)
            self.assertEqual(public.resolve(), app / "lcm")
            self.assertFalse(self.apps.exists())
        self.apps.mkdir()
        (self.apps / "study").mkdir()
        self.assertNotEqual(self.run_command("--application", "study", str(app), executable=self.installer).returncode, 0)
        self.assertTrue((self.apps / "study").is_dir())
        self.assertEqual(public.resolve(), app / "lcm")
        (self.apps / "study").rmdir()
        for _ in range(2):
            self.assertEqual(self.run_command("--application", "study", str(app), executable=self.installer).returncode, 0)
            self.assertEqual(public.resolve(), self.launcher)
            self.assertEqual((self.apps / "study").resolve(), app)

    def test_exec_preserves_pid_and_signal_delivery(self):
        target = self.app("study") / "lcm"
        target.write_text("#!/usr/bin/env python3\nimport os,signal,sys\n"
                          "signal.signal(signal.SIGTERM,lambda *_:sys.exit(29))\n"
                          "print(os.getpid(),flush=True)\nsignal.pause()\n")
        process = subprocess.Popen([str(self.launcher), "study"], cwd=self.cwd,
                                   env=self.env, stdout=subprocess.PIPE, text=True)
        try:
            with selectors.DefaultSelector() as ready:
                ready.register(process.stdout, selectors.EVENT_READ)
                self.assertTrue(ready.select(timeout=5), "application did not become ready")
                self.assertEqual(int(process.stdout.readline()), process.pid)
            process.send_signal(signal.SIGTERM)
            self.assertEqual(process.wait(timeout=5), 29)
        finally:
            if process.poll() is None:
                process.kill()
                process.wait(timeout=5)
            process.stdout.close()


if __name__ == "__main__":
    unittest.main()
