# Kubuntu workers: one setup, then start / stop

Already using the prepared `kubuntu` machine? Start with
[Run something now](../../TRY_IT.md). The steps below are for recreating its setup,
not a checklist you must repeat before each calculation.

The result is one private demo per Linux account:

```text
Browser → local gateway :8081 → private ts/SSH forwards → Kubuntu
                                                      ├─ TLS NATS broker
                                                      ├─ private HTTPS artifacts
                                                      ├─ demo-line: line calculations + separate REPLs
                                                      └─ demo-power: heavier calculations
```

Both worker profiles use isolated containers. Docker and Podman are alternative
local engines on the worker computer. Nothing installs Docker on your managed
browser computer. This recipe uses rootless engines, loopback listeners and the
single local `developer` identity; do not expose port 8081 publicly.
Anyone who can access that loopback endpoint acts as `developer`; this local mode
is not a login system or isolation between accounts sharing the gateway computer.

## 1. Connect to Kubuntu

**On the browser computer:**

```bash
ts ssh -tt kubuntu
```

You should get the Kubuntu shell. `-tt` gives you an interactive terminal for
any `sudo` password prompt. Passwords belong in that terminal, not in a config
file or a chat message. Later, `ts ssh -o BatchMode=yes kubuntu true` must work
without a password prompt so the demo can maintain its private tunnel.

## 2. Check the worker computer once

**On Kubuntu:**

```bash
podman info
systemctl --user is-system-running
```

Podman should report rootless operation and cgroup v2. A `degraded` user manager
can have an unrelated failed desktop unit; inspect failures rather than resetting
the whole account. Our two agents need an actual working user service manager.

On a **fresh Ubuntu/Kubuntu 24.04 host only**, install the prerequisites:

```bash
sudo apt-get update
sudo apt-get install podman crun uidmap slirp4netns fuse-overlayfs dbus-user-session openssl curl git
sudo loginctl enable-linger "$(id -un)"
```

This is the distribution-packaged Podman path documented by
[Podman](https://podman.io/docs/installation). Do not use `sudo podman` afterward.
Check that your account has non-overlapping subordinate UID/GID allocations:

```bash
rg "^$(id -un):" /etc/subuid /etc/subgid
```

If either is absent, have the administrator allocate a free range of at least
65,536 IDs; do not paste an arbitrary range over another account's allocation.

**Docker alternative:** install Docker Engine from its
[official Ubuntu instructions](https://docs.docker.com/engine/install/ubuntu/),
including `docker-ce-rootless-extras`, then follow its
[rootless setup](https://docs.docker.com/engine/security/rootless/):

```bash
dockerd-rootless-setuptool.sh install
docker context use rootless
systemctl --user start docker
docker info
```

Do not stop an existing rootful daemon if it serves other applications. On the
prepared Kubuntu machine rootless Docker is already installed, and the newly
installed rootful services were disabled. Choose **one engine for this demo**;
switching an existing instance is a deliberate reconfiguration, not a CLI toggle.

## 3. Put source and tools in one stable directory

**On Kubuntu**, the prepared layout is:

```text
/home/amauri/lcm-demo/
  source/                 this repository's source, without credentials
  tools/julia-1.12.7/      private Julia installation
  tools/caddy/caddy       private Caddy executable
  depot/                  Julia package cache
  state/                  generated private configuration and persistent data
```

For another account, replace `/home/amauri` throughout. Use paths without spaces
for this small operator helper. On a fresh host, clone your trusted repository
into `source` and check out the **same reviewed revision as the browser computer**.
Do not copy `.git` credentials, local `.env` files or runtime state from another
installation. Obtain Julia 1.12.7 from the
[official downloads and checksums](https://julialang.org/downloads/oldreleases/)
and a Caddy binary using the [official installation instructions](https://caddyserver.com/docs/install).
The prepared machine uses Caddy 2.10.2. Verify downloaded archives before extraction.

Set these variables in the Kubuntu shell:

```bash
export LCM_DEMO_ROOT=/home/amauri/lcm-demo
export JULIA_DEPOT_PATH="$LCM_DEMO_ROOT/depot"
export JULIA_LOAD_PATH=@:@stdlib
export PATH="$LCM_DEMO_ROOT/tools/julia-1.12.7/bin:$PATH"
cd "$LCM_DEMO_ROOT/source"
julia --startup-file=no --project=playground/runtime -e 'using Pkg; Pkg.instantiate()'
./playground/lcm runtime check-host --runtime podman
```

For Docker use `--runtime docker`. **Stop here if effective CPU, memory or PID
limits are missing.** The terminal must not fall back to an unrestricted native
REPL. Cgroup delegation is host administration; the
[Docker rootless resource-limits guide](https://docs.docker.com/engine/security/rootless/tips/#limiting-resources)
describes the systemd requirements. It is already configured on this Kubuntu host.

## 4. Install the three worker images and three infrastructure images

**On Kubuntu, from `source`:** skip building profiles that are already installed
and verified, as they are on the prepared machine. New builds can take time.

```bash
export LCM_DEMO_ENGINE=podman
for profile in julia-terminal line-parameters power-flow; do
    "$LCM_DEMO_ENGINE" build --target "$profile" \
        -f playground/worker/containers/Containerfile \
        -t "localhost/lcm-$profile:demo" .
    "$LCM_DEMO_ENGINE" image inspect "localhost/lcm-$profile:demo" --format '{{json .RepoDigests}}'
done
```

For Docker set `LCM_DEMO_ENGINE=docker`. Keep the returned **RepoDigest** for each
profile (`name@sha256:…`), not the mutable `:demo` tag or image ID. Digest values
can differ between the engines. A builder that does not return a local RepoDigest
needs an explicit local image export; do not substitute `.Id` to bypass validation.

Pull the pinned infrastructure images used by this recipe:

```bash
"$LCM_DEMO_ENGINE" pull docker.io/library/nats:2.11-alpine@sha256:e4bf19f15fd3218814a4e3c9e0064e1334bd8aa20d5984b9f1a0afd084f8cc00
"$LCM_DEMO_ENGINE" pull quay.io/minio/minio:RELEASE.2025-09-07T16-13-09Z@sha256:14cea493d9a34af32f524e538b8346cf79f3321eff8e708c1e2960462bd8936e
"$LCM_DEMO_ENGINE" pull quay.io/minio/mc:RELEASE.2025-08-13T08-35-41Z@sha256:a7fe349ef4bd8521fb8497f55c6042871b2ae640607cf99d9bede5e9bdf11727
```

## 5. Generate and install the Kubuntu configuration

**Still on Kubuntu:**

```bash
cp playground/deploy/demo/remote.example.toml "$LCM_DEMO_ROOT/remote.toml"
nano "$LCM_DEMO_ROOT/remote.toml"
```

Set the paths, choose `podman` or `docker`, and paste your three RepoDigests.
The `root` entry must name a **new** `state` directory. The helper refuses an
existing one, existing demo container names and a non-rootless engine. It never
silently replaces a deployment.

Check that ports **14222**, **14443** and **19000** and the `lcm-demo-*` /
`lcm-agent-demo-*` unit names are free before proceeding:

```bash
ss -lnt '( sport = :14222 or sport = :14443 or sport = :19000 )'
systemctl --user list-unit-files 'lcm-demo-*' 'lcm-agent-demo-*'
julia --startup-file=no --project=playground/runtime \
    playground/deploy/demo/configure.jl remote "$LCM_DEMO_ROOT/remote.toml"
systemctl --user link "$LCM_DEMO_ROOT"/state/units/*.service "$LCM_DEMO_ROOT"/state/units/*.target
systemctl --user daemon-reload
```

Expected: generated files, linked units, **no running calculations**. The service
files use the runtime's managed-agent renderer, including mandatory post-stop
executor cleanup. No password goes into Git. Broker and storage passwords are
random and scoped; the coordinator cannot write artifacts, and each worker can
only write/delete its own artifact prefix.

## 6. Connect the browser computer

**Exit the Kubuntu shell. On the browser computer:** use your existing built
checkout. For a fresh checkout, `./playground/bootstrap.sh` installs the publisher
requirements and builds the site; also instantiate `playground/runtime`.

```bash
cd /home/amartins/Documents/KUL/LineCableModels-playground
julia --startup-file=no --project=playground/runtime -e 'using Pkg; Pkg.instantiate()'
umask 077
mkdir -p /home/amartins/.local/state/lcm-demo
ts ssh -o BatchMode=yes kubuntu \
    'tar -C /home/amauri/lcm-demo/state -cf - coordinator' \
    | tar -xf - -C /home/amartins/.local/state/lcm-demo
cp playground/deploy/demo/client.example.toml /home/amartins/.local/state/lcm-demo/client.toml
nano /home/amartins/.local/state/lcm-demo/client.toml
```

This copies **only** the coordinator bundle over SSH: its client key, CA,
password, read-only storage credentials and profile references. Never copy the
entire remote `state` directory to the browser computer. Protect both locations;
the gateway is trusted to use these credentials.

Edit the local paths and `ts` adapter/host alias, then:

```bash
julia --startup-file=no --project=playground/runtime \
    playground/deploy/demo/configure.jl client /home/amartins/.local/state/lcm-demo/client.toml
systemctl --user link /home/amartins/.local/state/lcm-demo/units/*.service /home/amartins/.local/state/lcm-demo/units/*.target
systemctl --user daemon-reload
./playground/lcm demo provision
./playground/lcm demo start
```

Expected: private streams verified and a URL on **8081**. Ports **24222** and
**24443** are local SSH forwards, not public listeners. TLS certificate and host
verification remain enabled inside the tunnel. Kubuntu's MinIO HTTP backend is
loopback-only behind private Caddy HTTPS; this also avoids the desktop account's
MinIO certificate-watcher/inotify limitation. The MinIO console is not published.

If you choose a nondefault local state path, set `LCM_DEMO_STATE` to it before
each `lcm demo …` command. The example default needs no environment variables.

## 7. Approve the workers, then run something

Open **<http://127.0.0.1:8081/runtime/control>**. In worker registration:

1. Select provisioned identity **demo-line** and click **Enroll as pending**.
2. Select its registration, choose state **approved**, and apply the change.
3. Repeat for **demo-power**.
4. Wait for both to report **online**. Approval is persisted in the local SQLite
   database; you do not repeat it after a normal stop/start.

Follow [Run something now](../../TRY_IT.md) for the calculation and REPL clicks.
Worker registration is not model preparation: preparation stays explicit inside
each isolated application run.

## Daily use, shutdown, and logs

**On the browser computer:**

```bash
./playground/lcm demo start
./playground/lcm demo status
./playground/lcm demo stop
```

No services are enabled at boot. Stop closes the gateway first, lets managed
agents clean up their executors, then stops demo infrastructure and the tunnel.
It does not uninstall either engine, touch the publisher on 8080, prune containers,
or delete stored data. Restart creates new live sessions, not restored REPLs.

```bash
journalctl --user -u lcm-demo-tunnel -u lcm-demo-gateway -n 60 --no-pager
ts ssh kubuntu 'journalctl --user -u lcm-demo-broker -u lcm-demo-storage-tls -u lcm-agent-demo-line -u lcm-agent-demo-power -n 60 --no-pager'
```

Configuration, registrations and artifacts persist; Julia session variables and
executor scratch do not. This is **not yet persistent project editing or shared
user authentication**. Keep this distinction when evaluating the demonstration.

### Certificate expiry and backup

Certificates expire **30 days after setup**. Check on Kubuntu:

```bash
openssl x509 -in /home/amauri/lcm-demo/state/certs/server-cert.pem -noout -enddate
```

Renewal requires stopping the demo, issuing a new private development certificate
set, updating the agent/coordinator client certificates and CA together, copying
only the refreshed coordinator bundle locally, then starting again. Do not
disable TLS verification or rerun the generator over existing state to get around
expiry. For unattended use, replace this development PKI with managed certificates.

Back up both private state directories **while stopped** using an encrypted,
access-controlled backup. They contain credentials and the SQLite database;
they must not enter Git or a public file share. The `tools`, `depot`, and immutable
images can be rebuilt. Data volumes have no automatic deletion/retention policy
in this demo: monitor disk usage before treating it as a long-lived installation.

### Beyond the demo

The longer references are optional from here: [configuration](../../runtime/CONFIGURATION.md),
[physical isolation acceptance](../../runtime/PHYSICAL_ACCEPTANCE.md), and
[authenticated proxy deployment](../../runtime/PROXY.md). Acceptance harnesses
are destructive **to their own fixtures** and are not the way to start this demo.
