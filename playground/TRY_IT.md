# Run something now

This is the short route for the **already configured Kubuntu demo**. You do not
need to rebuild images, install Docker, or read the operator reference first.

## 1. Start it — on your usual computer

```bash
cd /home/amartins/Documents/KUL/LineCableModels-playground
./playground/lcm demo start
```

Open **<http://127.0.0.1:8081/workbenches/>**. Keep `127.0.0.1` in the address;
this demo's request-origin policy uses that exact host. The publisher on 8080 is
separate and is not stopped or replaced.

## 2. Run a real calculation

1. Find **Cable study** and click **Start isolated run**. Wait for the workbench.
2. Select **Line parameters** in its navigation.
3. Expand **Worker selection and preparation**. Choose the `line-parameters` profile. Automatic
   placement is fine; alternatively pin `demo-line`. Click **Assign worker**.
4. Click **Prepare executor**. Wait for ready. The first preparation is cold;
   assignment alone does not load the model or allow a calculation.
5. Collapse **Worker selection and preparation** to reveal the inputs and plot.
   Click **Run calculation**. A curve and a **Current result** indication should
   appear, together with the input fingerprint and execution provenance.
6. Change a frequency or separation input. The previous curve becomes
   **Outdated result**. Click **Run calculation** again to replace it.

The numerical process is on Kubuntu. The browser-facing application stays on
your computer. Keep the workbench tab open while using its session.
The first workbench launch can take roughly a minute; the line executor's first
preparation took about 29 seconds in the verified setup. Neither wait is a solver
running inside the web server.

For a second example, select **OHL / UGC case**, assign `power-flow` to
`demo-power`, then prepare and run. Cold preparation of this heavier profile can
take several minutes. There is no need to do it just to try the first example.

## 3. Try the Julia terminal

Select **Workers and preparation** in the same workbench. Under **Private Julia
terminal · optional**, assign the `julia-terminal` profile to `demo-line`.
Then select **Julia terminal** in the navigation and click **Connect**.
Wait for the Julia prompt and enter:

```julia
sum(1:100)
```

Expected answer: **5050**. This is a separate, limited container, not a REPL
inside the web server and not the scientific executor's workspace. Its variables
and scratch files are disposable. It is a private operator demo, not a public
multi-user terminal service.

## 4. Stop and come back later

```bash
./playground/lcm demo stop
./playground/lcm demo start
```

Stop closes live sessions and their owned executors. It keeps the private
configuration, worker registrations, runtime database and stored artifacts.
Start a **new isolated run** on your next visit; old session URLs are not projects
you can resume. There is no automatic start at login or boot.

## If something is not ready

```bash
./playground/lcm demo status
```

Worker inventory and diagnostics: <http://127.0.0.1:8081/runtime/control>.
Both `demo-line` and `demo-power` should be approved and online. If a start fails,
inspect the relevant log instead of retrying a preparation indefinitely:

```bash
journalctl --user -u lcm-demo-gateway -n 40 --no-pager
ts ssh kubuntu 'journalctl --user -u lcm-agent-demo-line -n 40 --no-pager'
```

To recreate this arrangement on a fresh worker computer, follow
[the ordered Kubuntu setup](deploy/demo/README.md). It separates commands run on
Kubuntu from commands run on the browser computer. The generated development
certificates expire after **30 days**; this is not an unattended public deployment.
