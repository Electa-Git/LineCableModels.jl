"""Manual warmed operational-overhead measurement, without a wall-clock test gate.

python3 test/gauntlet/progress_overhead.py [output.json]
Ten interleaved triples: off, publisher, publisher plus an external PTY watcher.
The Julia worker uses fixed seeded MC over 48 frequencies, at least six seconds
per run after calibration. Watcher startup/compilation precedes each timed run.
No real FEM/PSCAD calls, reports, checkpoints or user campaign directories.
"""
import fcntl
import json
import os
from pathlib import Path
import pty
import select
import signal
import statistics
import struct
import subprocess
import sys
import tempfile
import termios
import threading
import time

ROOT = Path(__file__).resolve().parents[2]
JULIA = ["julia", "--startup-file=no", "--project=gauntlet"]


def main():
    destination = Path(sys.argv[1] if len(sys.argv)>1 else "/tmp/lcm-progress-overhead.json")
    with tempfile.TemporaryDirectory(prefix="gauntlet overhead ") as directory:
        parent=Path(directory)
        log=open(parent/"worker.log", "w+")
        worker=subprocess.Popen(JULIA+[str(ROOT/"test/gauntlet/fixtures/progress_workload.jl")],
            cwd=ROOT, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=log, text=True,
            env={**os.environ,"JULIA_NUM_THREADS":"1","OPENBLAS_NUM_THREADS":"1"})
        def response(prefix):
            while True:
                line=worker.stdout.readline().strip()
                if line.startswith(prefix): return line.split("\t")[1:]
                if worker.poll() is not None:
                    log.seek(0)
                    raise RuntimeError(log.read())
        rows=[]
        started=time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime())
        try:
            trials,version,threads,blas=response("READY")
            print(f"Warmed {trials} scans × 48 frequencies; Julia {version}, threads={threads}, BLAS={blas}",flush=True)
            for repetition in range(10):
                order=["off","publisher","watcher"]
                order=order[repetition%3:]+order[:repetition%3]
                for arm in order:
                    root=parent/f"{repetition}-{arm}";root.mkdir()
                    watcher=None;stop=threading.Event();master=slave=None
                    if arm=="watcher":
                        master,slave=pty.openpty()
                        fcntl.ioctl(slave,termios.TIOCSWINSZ,struct.pack("HHHH",24,100,0,0))
                        def drain():
                            while not stop.is_set():
                                if select.select([master],[],[],0.1)[0]:
                                    try: os.read(master,65536)
                                    except OSError: return
                        reader=threading.Thread(target=drain,daemon=True);reader.start()
                        watcher=subprocess.Popen(JULIA+[str(ROOT/"test/gauntlet/fixtures/progress_process.jl"),
                            "watch",str(root),"measure"],cwd=ROOT,stdin=slave,stdout=slave,stderr=slave,
                            env={**os.environ,"TERM":"xterm-256color","NO_COLOR":"1","JULIA_NUM_THREADS":"1"},start_new_session=True)
                        deadline=time.monotonic()+120
                        while not (root/"watch.ready").exists():
                            assert watcher.poll() is None and time.monotonic()<deadline,"watcher startup failed"
                            time.sleep(0.05)
                        time.sleep(1.2) # Initial inventory render and compilation stay outside measurement.
                    try:
                        worker.stdin.write(f"{'off' if arm=='off' else 'auto'}\t{root}\n");worker.stdin.flush()
                        seconds,revisions,digest=response("RESULT")
                        row=dict(repetition=repetition+1,arm=arm,seconds=float(seconds),publications=int(revisions),digest=digest)
                        assert not rows or rows[0]["digest"]==digest,"scientific outputs changed across arms"
                        rows.append(row);print(json.dumps(row),flush=True)
                    finally:
                        if watcher:
                            watcher.send_signal(signal.SIGINT)
                            try: watcher.wait(timeout=10)
                            except subprocess.TimeoutExpired: watcher.terminate();watcher.wait(timeout=10)
                            stop.set();reader.join();os.close(master);os.close(slave)
            summary={}
            for arm in ("off","publisher","watcher"):
                values=[r["seconds"] for r in rows if r["arm"]==arm]
                paired=[100*(r["seconds"]/next(b["seconds"] for b in rows if b["arm"]=="off" and b["repetition"]==r["repetition"])-1)
                    for r in rows if r["arm"]==arm]
                q=statistics.quantiles(values,n=4,method="inclusive")
                pq=statistics.quantiles(paired,n=4,method="inclusive")
                summary[arm]=dict(paired_iqr_percent=pq[2]-pq[0],min_paired_percent=min(paired),max_paired_percent=max(paired),median_seconds=statistics.median(values),iqr_seconds=q[2]-q[0],
                    min_seconds=min(values),max_seconds=max(values),median_paired_percent=statistics.median(paired),
                    paired_percent=paired)
            report=dict(started_utc=started,kernel=os.uname().sysname,architecture=os.uname().machine,julia=version,threads=int(threads),blas_threads=int(blas),trials=int(trials),frequencies=48,
                repetitions=10,scope="Warmed operational owned MC call plus lifecycle observation/publication; excludes campaign reporting/checkpoints and watcher startup",
                target="Approximately <=1% median operational overhead; measurement only, no assertion",rows=rows,summary=summary)
            destination.write_text(json.dumps(report,indent=2)+"\n")
            print(json.dumps(summary,indent=2),flush=True)
        finally:
            worker.stdin.close()
            try: worker.wait(timeout=10)
            except subprocess.TimeoutExpired: worker.terminate();worker.wait(timeout=10)
            log.close()


if __name__=="__main__": main()
