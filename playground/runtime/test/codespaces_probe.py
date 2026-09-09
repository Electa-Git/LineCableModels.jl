#!/usr/bin/env python3
"""Codespaces-only Docker feasibility; not production integration acceptance."""
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import uuid

IMAGE = 'docker.io/library/julia:1.12-bookworm@sha256:709daad7eccdb0363b080df203601988cff8c6a9581162668b76804cca407888'
LABEL = 'org.linecablemodels.feasibility'

def command(args, seconds=30, environment=None):
    # Keep attach stdin open until the process ends. Inheriting the SSH command's
    # EOF loses the TTY stream with Docker start --attach --interactive.
    read_end, write_end = os.pipe()
    try:
        return subprocess.run(args, stdin=read_end, capture_output=True, text=True,
                              timeout=seconds, env=environment)
    finally:
        os.close(read_end)
        os.close(write_end)

def output(result):
    if result.returncode:
        raise RuntimeError(f'Exit {result.returncode}: {result.stdout[-1024:]} {result.stderr[-1024:]}')
    return result.stdout.strip()

def main():
    if os.environ.get('CODESPACES') != 'true':
        raise SystemExit('Refusing to run outside Codespaces.')
    directory = Path(__file__).parent
    source = (directory / 'ContainerIsolation.jl').read_text()
    nonce = str(uuid.uuid4())
    report = {'scope': 'docker-feasibility-only', 'full_integration': 'deferred',
              'source_sha256': hashlib.sha256(source.encode()).hexdigest(), 'image': IMAGE, 'checks': {}}
    owned_names = []
    baseline = None
    endpoint = json.loads(output(command(['docker', 'context', 'inspect', '--format', '{{json .Endpoints.docker.Host}}'])))
    if not re.fullmatch(r'unix:///[A-Za-z0-9/_.-]+', endpoint):
        raise SystemExit('Only an explicit local Unix-socket engine is permitted.')
    config = tempfile.mkdtemp(prefix='docker-config-', dir=directory)
    environment = {'PATH': os.environ['PATH'], 'DOCKER_CONFIG': config}

    def docker(*args, seconds=30):
        return command(['docker', '--host', endpoint, *args], seconds, environment)

    def inventory():
        return set(output(docker('ps', '-aq', '--no-trunc')).splitlines())

    def inspect(identity):
        item = json.loads(output(docker('inspect', identity)))[0]
        if item['Config'].get('Labels', {}).get(LABEL) != nonce:
            raise RuntimeError('Ownership mismatch; refusing operation.')
        return item

    def create(case, body, tty=False, cpu=True):
        name = 'lcm-feasibility-' + nonce + '-' + case
        owned_names.append(name)
        args = ['create', '--name', name, '--label', LABEL + '=' + nonce,
                '--pull=never', '--read-only', '--network=none', '--ipc=private', '--cgroupns=private',
                '--pid=', '--uts=', '--user=1000:1000', '--cap-drop=ALL', '--security-opt=no-new-privileges',
                '--log-driver=none', '--restart=no', '--workdir=/tmp', '--interactive',
                '--memory=536870912', '--memory-swap=536870912', '--pids-limit=64', '--shm-size=65536',
                '--ulimit=core=0:0', '--ulimit=msgqueue=0:0', '--ulimit=nofile=1024:1024',
                '--tmpfs', '/tmp:rw,nosuid,nodev,noexec,size=8323072,mode=1777',
                '--entrypoint=/usr/local/julia/bin/julia']
        if cpu:
            args += ['--cpu-period=100000', '--cpu-quota=50000']
        if tty:
            args += ['--tty']
        for key, value in {'HOME': '/tmp', 'TMPDIR': '/tmp', 'JULIA_DEPOT_PATH': '/tmp/depot',
                           'JULIA_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'JULIA_LOAD_PATH': '@stdlib'}.items():
            args += ['--env', key + '=' + value]
        args += [IMAGE, '--startup-file=no', '--history-file=no', '-e', source + '\n' + body]
        identity = output(docker(*args))
        if not re.fullmatch('[a-f0-9]{64}', identity):
            raise RuntimeError('Invalid container ID.')
        item = inspect(identity)
        host = item['HostConfig']
        if host.get('Binds') or host.get('Devices') or host.get('Privileged'):
            raise RuntimeError('Unexpected mounts, devices or privilege.')
        allowed = {'PATH', 'JULIA_PATH', 'JULIA_VERSION', 'JULIA_GPG', 'HOME', 'TMPDIR',
                   'JULIA_DEPOT_PATH', 'JULIA_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'JULIA_LOAD_PATH'}
        if any(entry.split('=', 1)[0] not in allowed for entry in item['Config']['Env']):
            raise RuntimeError('Unexpected environment; refusing start.')
        report['checks'][case + '_policy'] = {key: host.get(key) for key in
            ('ReadonlyRootfs', 'NetworkMode', 'CgroupnsMode', 'Memory', 'MemorySwap', 'PidsLimit',
             'CpuQuota', 'CpuPeriod', 'CapDrop', 'SecurityOpt')}
        return identity

    try:
        version = json.loads(output(docker('version', '--format', '{{json .Server}}')))
        if 'podman' in json.dumps(version).lower():
            raise RuntimeError('Podman shim is not genuine Docker evidence.')
        info = json.loads(output(docker('info', '--format', '{{json .}}')))
        report['engine'] = {key: info.get(key) for key in ('ServerVersion', 'OSType', 'CgroupVersion',
            'CgroupDriver', 'MemoryLimit', 'SwapLimit', 'CpuCfsPeriod', 'CpuCfsQuota', 'PidsLimit', 'SecurityOptions')}
        if info.get('OSType') != 'linux' or info.get('CgroupVersion') != '2' or not all(
                info.get(key) is True for key in ('MemoryLimit', 'SwapLimit', 'CpuCfsPeriod', 'CpuCfsQuota', 'PidsLimit')):
            raise RuntimeError('Host resource prerequisites unavailable.')
        report['pid1'] = Path('/proc/1/comm').read_text().strip()
        manager = command(['busctl', '--user', '--timeout=5', 'get-property', 'org.freedesktop.systemd1',
                           '/org/freedesktop/systemd1', 'org.freedesktop.systemd1.Manager', 'Version'])
        report['user_systemd_manager'] = manager.returncode == 0
        report['managed_agent_acceptance'] = 'not-tested' if manager.returncode == 0 else 'unavailable'
        baseline = inventory()
        image = json.loads(output(docker('image', 'inspect', IMAGE)))[0]
        report['image_id'] = image['Id']
        report['repository_digests'] = image['RepoDigests']
        guard = 'limits = ContainerLimits(0.5, 536870912, 64, 8388608)\n'
        body = guard + '''
println("ISOLATION_OK ", verify_container_isolation(limits))
println("JULIA_VERSION ", VERSION)
println("JULIA_RESULT ", sum(1:100))
write("/tmp/feasibility.txt", "disposable scratch")
@assert read("/tmp/feasibility.txt", String) == "disposable scratch"
root_denied = try
    write("/lcm-forbidden-write", "forbidden"); false
catch
    true
end
@assert root_denied
@assert !ispath("/var/run/docker.sock")
before = read("/sys/fs/cgroup/cpu.stat", String)
deadline = time() + 1.5
while time() < deadline
    sin(time())
end
after = read("/sys/fs/cgroup/cpu.stat", String)
throttles(s) = parse(Int, match(r"nr_throttled (\\d+)", s).captures[1])
@assert throttles(after) > throttles(before)
println("CPU_THROTTLING_OBSERVED")
println("SCRATCH_OK ROOT_WRITE_DENIED NO_ENGINE_SOCKET")
'''
        for case, tty in (('plain', False), ('tty', True)):
            identity = create(case, body, tty=tty)
            run = docker('start', '--attach', '--interactive', identity, seconds=120)
            report['checks'][case] = {'status': 'incomplete', 'exit_code': run.returncode,
                                    'output': run.stdout, 'stderr': run.stderr,
                                    'container_state': inspect(identity)['State']}
            output(run)
            if any(marker not in run.stdout for marker in ('ISOLATION_OK', 'JULIA_RESULT 5050', 'CPU_THROTTLING_OBSERVED')):
                raise RuntimeError('Incomplete Julia evidence.')
            report['checks'][case]['status'] = 'pass'
            print(case + ': Julia and production isolation verifier passed.', flush=True)
        negative = guard + '''
try
    verify_container_isolation(limits)
catch failure
    failure isa IsolationError && failure.code == :cpu_limit_unverified || rethrow()
    println("EXPECTED_REJECTION cpu_limit_unverified")
    exit(78)
end
println("UNSAFE_USER_CODE_REACHED")
'''
        identity = create('missing-cpu', negative, cpu=False)
        run = docker('start', '--attach', '--interactive', identity, seconds=120)
        state = inspect(identity)['State']
        if state['ExitCode'] != 78 or 'EXPECTED_REJECTION cpu_limit_unverified' not in run.stdout or 'UNSAFE_USER_CODE_REACHED' in run.stdout:
            raise RuntimeError('Missing CPU limit did not fail closed.')
        report['checks']['missing_cpu'] = {'status': 'pass', 'exit_code': state['ExitCode'], 'output': run.stdout.strip()}
        identity = create('cancel', guard + 'verify_container_isolation(limits); sleep(60)')
        output(docker('start', identity))
        inspect(identity)
        output(docker('kill', '--signal=KILL', identity))
        if inspect(identity)['State']['Running']:
            raise RuntimeError('Canceled container still running.')
        report['checks']['cancel'] = {'status': 'pass'}
        report['result'] = 'docker-feasibility-pass'
    except Exception as failure:
        report['result'] = 'failed'
        report['error'] = str(failure)[:4096]
    finally:
        errors = []
        for name in owned_names:
            try:
                if output(docker('ps', '-aq', '--filter', 'name=^/' + name + '$')):
                    inspect(name)
                    output(docker('rm', '--force', name))
            except Exception as failure:
                errors.append(str(failure)[:1024])
        report['cleanup_errors'] = errors
        report['baseline_restored'] = baseline is not None and inventory() == baseline
        if errors or not report['baseline_restored']:
            report['result'] = 'failed'
        (directory / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps(report, indent=2), flush=True)
    return 0 if report['result'] == 'docker-feasibility-pass' else 1

if __name__ == '__main__':
    raise SystemExit(main())
