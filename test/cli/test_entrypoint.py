"""Exercise the package launcher without Julia packages or application services."""
import os
from pathlib import Path
import shutil
import signal
import subprocess
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[2]


class EntryPointTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="lcm cli ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.root = self.directory / "package checkout"
        self.root.mkdir()
        shutil.copytree(REPOSITORY / "cli", self.root / "cli")
        (self.root / "Project.toml").write_text('name = "LineCableModels"\nversion = "0.2.0"\n')
        self.bin = self.directory / "user bin"
        self.bin.mkdir()
        self.links = self.directory / "configuration" / "applications"
        self.env = dict(os.environ, LCM_APPLICATIONS_DIR=str(self.links),
                        LCM_INSTALL_DIR=str(self.bin), PATH=f"{self.bin}:{os.environ['PATH']}")
        self.caller = self.directory / "working directory"
        self.caller.mkdir()
        self.launcher = self.root / "cli" / "lcm"
        self.installer = self.root / "cli" / "install.sh"

    def invoke(self, *arguments, launcher=None, **kwargs):
        return subprocess.run([str(launcher or self.launcher), *arguments],
                              cwd=self.caller, env=self.env, capture_output=True, **kwargs)

    def application(self, name, directory=None):
        target = (directory or self.root) / name
        target.mkdir(parents=True)
        script = target / "lcm"
        script.write_text('#!/usr/bin/env bash\n'
                          'printf "%s\\0" "$PWD" "$0" "$@"\n'
                          'exit "${LCM_TEST_EXIT:-0}"\n')
        script.chmod(0o755)
        return target

    def test_global_inspection_does_not_start_an_application(self):
        self.application("gauntlet")
        self.application("playground")
        for arguments in [(), ("--help",), ("-h",), ("help",), ("--paths",)]:
            with self.subTest(arguments=arguments):
                result = self.invoke(*arguments)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn(b"gauntlet", result.stdout)
                self.assertIn(b"playground", result.stdout)
                self.assertNotIn(b"\0", result.stdout)
                self.assertNotIn(b"  cli ", result.stdout)
        self.assertEqual(self.invoke("--version").stdout, b"lcm 0.2.0\n")

    def test_symlink_chain_arguments_cwd_and_application_exit_status(self):
        app = self.application("gauntlet")
        intermediate = self.bin / "intermediate"
        intermediate.symlink_to(self.launcher)
        public = self.bin / "lcm"
        public.symlink_to("intermediate")
        self.env["LCM_TEST_EXIT"] = "37"
        arguments = ("gauntlet", "run", "--definition", "file with spaces.jl", "", "$(touch forbidden)")
        result = self.invoke(*arguments, launcher=public)
        self.assertEqual(result.returncode, 37)
        self.assertEqual(result.stdout.decode().split("\0")[:-1],
                         [str(self.caller), str(app / "lcm"), *arguments])
        self.assertFalse((self.caller / "forbidden").exists())
        result = self.invoke("--paths", launcher="lcm")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(str(self.root).encode(), result.stdout)

    def test_unknown_invalid_and_global_extra_arguments(self):
        for arguments in [("unknown",), ("../cli",), ("cli",), ("--bogus",),
                          ("",), ("--version", "unexpected"), ("help", "a", "b")]:
            with self.subTest(arguments=arguments):
                result = self.invoke(*arguments)
                self.assertEqual(result.returncode, 2)
                self.assertTrue(result.stderr)

    def test_other_worktree_and_playground_commands(self):
        app = self.application("playground", self.directory / "other worktree")
        result = self.invoke("--application", "playground", str(app), launcher=self.installer)
        self.assertEqual(result.returncode, 0, result.stderr)
        for command in ("playground", "runtime", "worker", "presentation", "nats", "container", "demo"):
            with self.subTest(command=command):
                result = self.invoke(command, "--help", launcher=self.bin / "lcm")
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(result.stdout.decode().split("\0")[-3:-1], [command, "--help"])
        self.assertIn(str(app).encode(), self.invoke("--paths").stdout)
        self.assertEqual(self.invoke("help", "playground").stdout,
                         self.invoke("playground", "--help").stdout)

    def test_checkout_application_wins_and_future_application_uses_same_contract(self):
        configured = self.application("study", self.directory / "another checkout")
        self.links.mkdir(parents=True)
        (self.links / "study").symlink_to(configured, target_is_directory=True)
        local = self.application("study")
        result = self.invoke("study", "inspect")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(str(local / "lcm").encode(), result.stdout)
        paths = self.invoke("--paths").stdout
        self.assertIn(str(local).encode(), paths)
        self.assertNotIn(str(configured).encode(), paths)

    def test_gauntlet_owns_julia_project_and_consumes_its_namespace(self):
        app = self.root / "gauntlet"
        app.mkdir()
        shutil.copy2(REPOSITORY / "gauntlet" / "lcm", app / "lcm")
        julia = self.application("fake julia", self.directory) / "lcm"
        self.env["LCM_JULIA"] = str(julia)
        result = self.invoke("gauntlet", "run", "--definition", "relative case.jl")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.decode().split("\0")[:-1], [str(self.caller), str(julia),
                         f"--project={app}", "--startup-file=no", str(app / "cli.jl"),
                         "run", "--definition", "relative case.jl"])

    def test_installer_is_repeatable_and_replaces_legacy_symlink(self):
        old = self.application("playground", self.directory / "legacy checkout")
        (self.bin / "lcm").symlink_to(old / "lcm")
        for _ in range(2):
            result = self.invoke("--application", "playground", str(old), launcher=self.installer)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((self.bin / "lcm").resolve(), self.launcher)
            self.assertEqual((self.links / "playground").resolve(), old)

    def test_installer_validates_all_inputs_before_changing_links(self):
        existing = self.application("gauntlet", self.directory / "existing")
        public = self.bin / "lcm"
        public.symlink_to(existing / "lcm")
        for arguments in [("--application", "missing", str(self.directory)),
                          ("--application", "../bad", str(existing)),
                          ("--application", "gauntlet"), ("--bin-dir",),
                          ("--application", "runtime", str(existing)),
                          ("--application", "gauntlet", str(existing),
                           "--application", "gauntlet", str(existing))]:
            with self.subTest(arguments=arguments):
                result = self.invoke(*arguments, launcher=self.installer)
                self.assertNotEqual(result.returncode, 0)
                self.assertEqual(public.resolve(), existing / "lcm")
                self.assertFalse(self.links.exists())

    def test_installer_preserves_regular_files_and_directories(self):
        public = self.bin / "lcm"
        public.write_text("owned file")
        self.assertNotEqual(self.invoke(launcher=self.installer).returncode, 0)
        self.assertEqual(public.read_text(), "owned file")
        public.unlink()
        self.links.mkdir(parents=True)
        (self.links / "gauntlet").mkdir()
        app = self.application("gauntlet")
        self.assertNotEqual(self.invoke("--application", "gauntlet", str(app),
                                        launcher=self.installer).returncode, 0)
        self.assertTrue((self.links / "gauntlet").is_dir())
        self.assertFalse(public.exists())

    def test_application_receives_termination_directly(self):
        app = self.application("study")
        (app / "lcm").write_text('#!/usr/bin/env bash\n'
                                 "trap 'exit 42' TERM\n"
                                 "printf 'ready\\n'\n"
                                 'read -r ignored\n')
        process = subprocess.Popen([str(self.launcher), "study"], cwd=self.caller,
                                   env=self.env, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        try:
            self.assertEqual(process.stdout.readline(), b"ready\n")
            process.send_signal(signal.SIGTERM)
            self.assertEqual(process.wait(timeout=5), 42)
        finally:
            if process.poll() is None:
                process.kill()
                process.wait()
            process.stdin.close()
            process.stdout.close()


if __name__ == "__main__":
    unittest.main()
