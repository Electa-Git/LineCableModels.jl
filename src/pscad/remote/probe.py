import os
import sys
import threading
import traceback


def timed_out():
    try:
        print(
            f"PSCAD remote-call probe exceeded {timeout_seconds:g} s",
            file=sys.stderr,
            flush=True,
        )
    finally:
        os._exit(124)


try:
    host, port_text, timeout_text = sys.argv[1:]
    port = int(port_text)
    timeout_seconds = float(timeout_text)

    # A watchdog bounds imports and remote calls. Each exit bypasses interpreter shutdown.
    watchdog = threading.Timer(timeout_seconds, timed_out)
    watchdog.daemon = True
    watchdog.start()

    print(f"Probing PSCAD on {host}:{port} (limit={timeout_seconds:g} s)", flush=True)
    import mhi.pscad

    app = mhi.pscad.connect(host=host, port=port, timeout=min(30, timeout_seconds))
    busy = app.is_busy()
    version = app.version
    print(f"PSCAD remote calls answered: busy={busy}, version={version}", flush=True)
except BaseException:
    try:
        traceback.print_exc(file=sys.stderr)
        sys.stderr.flush()
    finally:
        os._exit(1)
os._exit(0)
