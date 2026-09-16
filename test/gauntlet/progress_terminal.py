"""Manual PTY screen smoke test; runs two Julia processes and no native solvers.

Run from the repository root: python3 test/gauntlet/progress_terminal.py
Uses only the Python standard library. The small VT emulator checks the rendered
screen, including cursor updates; escape-sequence presence alone is insufficient.
"""
import codecs
import fcntl
import os
from pathlib import Path
import pty
import re
import select
import signal
import struct
import subprocess
import tempfile
import termios
import time
import tomllib
import unicodedata

ROOT = Path(__file__).resolve().parents[2]
FIXTURE = ROOT / "test/gauntlet/fixtures/progress_process.jl"


class Screen:
    def __init__(self, height=24, width=100):
        self.height, self.width = height, width
        self.rows = [[" "] * width for _ in range(height)]
        self.row = self.col = 0
        self.pending = ""
        self.visible = True
        self.decoder = codecs.getincrementaldecoder("utf-8")("replace")

    def resize(self, height, width):
        # Hard line endings are retained; only soft terminal wraps can reflow.
        self.rows = [(r + [" "] * width)[:width] for r in self.rows[:height]]
        self.rows += [[" "] * width for _ in range(height - len(self.rows))]
        self.height, self.width = height, width
        self.row, self.col = min(self.row, height - 1), min(self.col, width - 1)

    def feed(self, data):
        self.pending += self.decoder.decode(data)
        while self.pending:
            if self.pending.startswith("\x1b]"):
                end=self.pending.find("\x07")
                if end<0: break
                self.pending=self.pending[end+1:] # Julia sets the terminal title at startup.
                continue
            if self.pending.startswith("\x1b"):
                match = re.match(r"\x1b\[([?0-9;]*)([ -/]*)([@-~])", self.pending)
                if not match:
                    break
                self.pending = self.pending[match.end():]
                params, _, code = match.groups()
                count = int(params or "1") if params.isdigit() or not params else 1
                if code == "A": self.row = max(0, self.row - count)
                elif code == "B": self.row = min(self.height - 1, self.row + count)
                elif code == "C": self.col = min(self.width - 1, self.col + count)
                elif code == "G": self.col = min(self.width - 1, max(0, count - 1))
                elif code == "K":
                    assert params == "2", match.group()
                    self.rows[self.row] = [" "] * self.width
                elif code == "J" and params in ("", "0"):
                    self.rows[self.row][self.col:] = [" "] * (self.width-self.col)
                    for row in range(self.row+1,self.height):
                        self.rows[row] = [" "] * self.width
                elif code == "m": pass
                elif code in "hl" and params == "?25": self.visible = code == "h"
                else: raise AssertionError(f"Unhandled VT sequence {match.group()!r}")
                continue
            char, self.pending = self.pending[0], self.pending[1:]
            if char == "\r": self.col = 0
            elif char == "\n":
                self.row += 1
                if self.row == self.height:
                    self.rows.pop(0)
                    self.rows.append([" "] * self.width)
                    self.row -= 1
            elif char >= " ":
                width = 0 if unicodedata.combining(char) else (2 if unicodedata.east_asian_width(char) in "WF" else 1)
                if width:
                    if self.col >= self.width:
                        self.feed(b"\r\n")
                    self.rows[self.row][self.col] = char
                    self.col += width

    def text(self):
        return ["".join(row).rstrip() for row in self.rows if any(c != " " for c in row)]


def run_smoke():
    with tempfile.TemporaryDirectory(prefix="gauntlet terminal '") as directory:
        trace = Path(tempfile.mkdtemp(prefix="lcm-progress-provisional-")) / "terminal.raw"
        print(f"PTY evidence: {trace}", flush=True)
        master, slave = pty.openpty()
        screen = Screen()
        def resize(height, width):
            fcntl.ioctl(slave, termios.TIOCSWINSZ, struct.pack("HHHH", height, width, 0, 0))
            screen.resize(height, width)
        resize(24, 100)
        env = {**os.environ, "TERM": "xterm-256color", "NO_COLOR": "1"}
        command = ["julia", "--startup-file=no", "--compiled-modules=existing", "--project=gauntlet", str(FIXTURE)]
        producer = subprocess.Popen(command + ["produce", directory], cwd=ROOT, env=env,
                                    stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        watcher = subprocess.Popen(command + ["watch", directory, "current-session"], cwd=ROOT, env=env,
                                   stdin=slave, stdout=slave, stderr=slave, start_new_session=True)
        def pump(seconds):
            deadline = time.monotonic() + seconds
            while time.monotonic() < deadline:
                if select.select([master], [], [], min(0.1, max(0, deadline-time.monotonic())))[0]:
                    data=os.read(master, 65536)
                    with trace.open("ab") as stream: stream.write(data)
                    screen.feed(data)
        def until(predicate, timeout=60):
            deadline = time.monotonic() + timeout
            while not predicate():
                assert time.monotonic() < deadline, (screen.text(), screen.pending[:500])
                assert watcher.poll() is None, screen.text()
                assert producer.poll() in (None, 0), producer.stderr.read().decode()
                pump(0.1)
        def command_phase(name): (Path(directory) / name).touch()
        def text(): return "\n".join(screen.text())
        try:
            until(lambda: "expanded current label" in text())
            assert len(screen.text()) == 6, screen.text()
            command_phase("compact")
            until(lambda: "compact" in text() and "expanded current label" not in text())
            assert len(screen.text()) == 6 and "wide label" not in text(), screen.text()
            resize(24, 45); pump(1.5)
            assert len(screen.text()) <= 1, screen.text()
            resize(24, 100); pump(1.5)
            assert len(screen.text()) == 6, screen.text()
            command_phase("busy")
            until(lambda: "solving" in text())
            first = screen.text()[4]
            pump(0.3)
            assert screen.text()[4] != first, screen.text()
            pump(4.0)
            command_phase("sample")
            until(lambda: "observation suspended" in text())
            assert "scans" not in screen.text()[4], screen.text()
            before = screen.text()[4]
            pump(0.3)
            assert screen.text()[4] != before, screen.text()
            until(lambda: "observation suspended" not in text())
            command_phase("finish")
            until(lambda: "ETA done" in text())
            assert "Complete 1" in text() and "Failed 0" in text() and len(screen.text()) == 6, screen.text()
            final = screen.text()
            pump(1.5)
            assert screen.text() == final, screen.text()
            assert producer.wait(timeout=10) == 0, producer.stderr.read().decode()
            watcher.send_signal(signal.SIGINT)
            pump(1.0)
            assert watcher.wait(timeout=10) == 0
            assert screen.visible
            print("PTY screen smoke passed: six rows, shorter text, narrow resize, opaque spinner, suspended sample, frozen closure, cursor restored.")
            print("\n".join(final))
        finally:
            for process in (watcher, producer):
                if process.poll() is None:
                    process.terminate()
                    process.wait(timeout=10)
            os.close(master); os.close(slave)


if __name__ == "__main__":
    run_smoke()
