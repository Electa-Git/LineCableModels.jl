using Test, LineCableModelsExecutionCore
const IsolationCore = LineCableModelsExecutionCore

struct IsolationStat
    uid::Int
    mode::UInt
    rdev::Int
    kind::Symbol
end

function test_rootless_device_mounts()
@testset "rootless OCI character-device binds are exact, not arbitrary mounts" begin
    devices = ("/dev/null"=>0x0103, "/dev/zero"=>0x0105, "/dev/full"=>0x0107,
        "/dev/random"=>0x0108, "/dev/urandom"=>0x0109, "/dev/tty"=>0x0500)
    for (path, number) in devices, uid in (0,65534)
        f = isolation_fixture()
        f.files["/proc/self/mountinfo"] *= "\n" * f.mount(path,"devtmpfs","rw,nosuid")
        f.stats[path] = IsolationStat(uid,0o666,number,:device)
        @test f.check().uid == 1000
        for (owner, device, kind) in ((uid,number,:file),(uid,number+1,:device),(1000,number,:device))
            f.stats[path] = IsolationStat(owner,0o666,device,kind)
            @test_throws IsolationError f.check()
        end
    end
    f = isolation_fixture()
    f.files["/proc/self/mountinfo"] = replace(f.files["/proc/self/mountinfo"],
        "/proc/kcore rw,nosuid - tmpfs"=>"/proc/kcore rw,nosuid - devtmpfs")
    f.stats["/proc/kcore"] = IsolationStat(65534,0o666,0x0103,:device)
    f.stats["/dev/null"] = IsolationStat(65534,0o666,0x0103,:device)
    @test f.check().uid == 1000
    for number in (0x0101,0x0105,0x0800)
        f.stats["/proc/kcore"] = IsolationStat(65534,0o666,number,:device)
        f.stats["/dev/null"] = IsolationStat(65534,0o666,number,:device)
        @test_throws IsolationError f.check()
    end
    for fs in ("ext4","tmpfs"), options in ("rw,nosuid","rw")
        f = isolation_fixture()
        f.files["/proc/self/mountinfo"] *= "\n" * f.mount("/dev/zero",fs,options)
        f.stats["/dev/zero"] = IsolationStat(65534,0o666,0x0105,:device)
        @test_throws IsolationError f.check()
    end
end
end
Base.isdir(s::IsolationStat) = s.kind == :directory
Base.isfile(s::IsolationStat) = s.kind == :file
Base.ischardev(s::IsolationStat) = s.kind == :device

function isolation_fixture()
    status = "Uid:\t1000 1000 1000 1000\nGid:\t1000 1000 1000 1000\nGroups:\t1000\n" *
        join((key * ":\t0000000000000000\n" for key in ("CapInh", "CapPrm", "CapEff", "CapBnd", "CapAmb"))) *
        "NoNewPrivs:\t1\nSeccomp:\t2\n"
    mount(path, fs, options) = "1 0 0:1 / $path $options - $fs none $options"
    files = Dict(
        "/proc/self/status"=>status, "/proc/1/status"=>status,
        "/proc/self/cgroup"=>"0::/\n", "/sys/fs/cgroup/cpu.max"=>"100000 100000\n",
        "/sys/fs/cgroup/memory.max"=>"1073741824\n", "/sys/fs/cgroup/memory.swap.max"=>"0\n",
        "/sys/fs/cgroup/pids.max"=>"128\n", "/proc/net/dev"=>"Inter-| Receive | Transmit\n lo: 0 0 0 0\n",
        "/proc/self/limits"=>"Max open files            1024 1024 files\nMax msgqueue size          0 0 bytes\nMax core file size         0 0 bytes\n",
        "/proc/self/mountinfo"=>join([
            mount("/", "overlay", "ro"), mount("/proc", "proc", "rw,nosuid,nodev,noexec"),
            mount("/sys", "sysfs", "ro,nosuid,nodev,noexec"), mount("/sys/fs/cgroup", "cgroup2", "ro,nosuid,nodev,noexec"),
            mount("/tmp", "tmpfs", "rw,nosuid,nodev,noexec"), mount("/dev", "tmpfs", "rw,nosuid"),
            mount("/dev/pts", "devpts", "rw,nosuid,noexec"), mount("/dev/shm", "tmpfs", "rw,nosuid,nodev,noexec"),
            mount("/dev/mqueue", "mqueue", "rw,nosuid,nodev,noexec"),
            mount("/etc/hosts", "ext4", "rw"), mount("/etc/hostname", "ext4", "rw"), mount("/etc/resolv.conf", "ext4", "rw"),
            mount("/proc/kcore", "tmpfs", "rw,nosuid"), mount("/proc/sys", "proc", "ro,nosuid,nodev,noexec"),
        ], '\n'))
    stats = Dict("/dev"=>IsolationStat(0,0o755,0,:directory),
        "/dev/null"=>IsolationStat(0,0o666,259,:device), "/proc/kcore"=>IsolationStat(0,0o666,259,:device))
    for path in ("/etc/hosts", "/etc/hostname", "/etc/resolv.conf")
        stats[path] = IsolationStat(0,0o644,0,:file)
    end
    disks = Dict("/tmp"=>(ftype=0x01021994,bsize=4096,blocks=65520),
        "/dev/shm"=>(ftype=0x01021994,bsize=4096,blocks=16))
    limits = ContainerLimits(1,1024^3,128,256*1024^2)
    check() = IsolationCore.check_container_isolation(limits, p->files[p], p->stats[p], p->disks[p])
    return (; files, stats, disks, limits, check, mount)
end

@testset "effective container policy: bounded kernel evidence" begin
    f = isolation_fixture()
    @test f.check() == (schema_version=1,cpus=1.0,memory_bytes=1024^3,pids=128,scratch_bytes=256*1024^2,uid=1000)
    for args in ((0.001,1024^3,128,1024^2),(true,1024^3,128,1024^2),
            (Inf,1024^3,128,1024^2),(1,true,128,1024^2),(1,1024^3,65537,1024^2),
            (1,1024^3,128,65536))
        @test_throws ArgumentError ContainerLimits(args...)
    end
    for (path, value) in (
            "/proc/self/cgroup"=>"0::/host/other", "/sys/fs/cgroup/cpu.max"=>"max 100000",
            "/sys/fs/cgroup/cpu.max"=>"100001 100000", "/sys/fs/cgroup/cpu.max"=>"1 1",
            "/sys/fs/cgroup/memory.max"=>"1073741825", "/sys/fs/cgroup/memory.max"=>"max",
            "/sys/fs/cgroup/memory.swap.max"=>"1", "/sys/fs/cgroup/pids.max"=>"129",
            "/sys/fs/cgroup/pids.max"=>"max", "/proc/net/dev"=>"lo: 0\neth0: 0\n",
            "/proc/net/dev"=>"", "/proc/self/mountinfo"=>repeat("x",256*1024+1))
        f = isolation_fixture(); f.files[path] = value
        @test_throws IsolationError f.check()
    end
    for path in ("/proc/self/status", "/proc/1/status"), (old,new) in (
            "Uid:\t1000 1000 1000 1000"=>"Uid:\t0 0 0 0", "Groups:\t1000"=>"Groups:\t1000 27",
            "CapBnd:\t0000000000000000"=>"CapBnd:\t0000000000000001",
            "NoNewPrivs:\t1"=>"NoNewPrivs:\t0", "Seccomp:\t2"=>"Seccomp:\t0")
        f = isolation_fixture(); f.files[path] = replace(f.files[path],old=>new)
        @test_throws IsolationError f.check()
    end
    for (old,new) in ("1024 1024"=>"1024 unlimited", "1024 1024"=>"1025 1025",
            "0 0 bytes"=>"0 1 bytes")
        f = isolation_fixture(); f.files["/proc/self/limits"] = replace(f.files["/proc/self/limits"],old=>new)
        @test_throws IsolationError f.check()
    end
    for (path,fs,options) in (("/owned", "tmpfs", "rw"), ("/run/secrets", "tmpfs", "ro"),
            ("/dev/socket", "ext4", "rw"), ("/tmp", "tmpfs", "rw,nosuid,nodev,noexec"))
        f = isolation_fixture(); f.files["/proc/self/mountinfo"] *= "\n" * f.mount(path,fs,options)
        @test_throws IsolationError f.check()
    end
    for path in ("/etc/hosts", "/dev", "/proc/kcore")
        f = isolation_fixture(); f.stats[path] = IsolationStat(1000,0o666,259,:device)
        @test_throws IsolationError f.check()
    end
    for (old,new) in (("/ / ro - overlay none ro", "/ / rw - overlay none rw"), ("/ /tmp rw,nosuid,nodev,noexec", "/ /tmp rw,nosuid,nodev"))
        f = isolation_fixture(); f.files["/proc/self/mountinfo"] = replace(f.files["/proc/self/mountinfo"],old=>new)
        @test_throws IsolationError f.check()
    end
    f = isolation_fixture(); f.disks["/tmp"] = (ftype=0x01021994,bsize=4096,blocks=65536)
    @test_throws IsolationError f.check()
    f = isolation_fixture(); f.disks["/dev/shm"] = (ftype=0x01021994,bsize=4096,blocks=17)
    @test_throws IsolationError f.check()
    f = isolation_fixture(); f.disks["/tmp"] = (ftype=0x1234,bsize=4096,blocks=1)
    @test_throws IsolationError f.check()
    f = isolation_fixture(); f.stats["/proc/kcore"] = IsolationStat(0,0o666,260,:device)
    @test_throws IsolationError f.check()
    f = isolation_fixture(); f.files["/proc/self/status"] = replace(f.files["/proc/self/status"],"Groups:\t1000\n"=>"")
    @test_throws IsolationError f.check()
    f = isolation_fixture(); f.files["/proc/self/status"] *= "Uid:\t1000 1000 1000 1000\n"
    @test_throws IsolationError f.check()
    # A terminal's devpts-backed console is an expected mount, not scratch storage.
    f = isolation_fixture(); f.stats["/dev/console"] = IsolationStat(1000,0o620,3,:device)
    f.files["/proc/self/mountinfo"] *= "\n" * f.mount("/dev/console","devpts","rw,nosuid,noexec")
    @test f.check().uid == 1000
    @test !occursin("/proc", sprint(showerror,IsolationError(:network_not_isolated)))
end

test_rootless_device_mounts()

@testset "container identity metadata is read-only and root-owned" begin
    for path in ("/etc/passwd","/etc/group","/run/.containerenv")
        f = isolation_fixture()
        f.files["/proc/self/mountinfo"] *= "\n" * f.mount(path,"tmpfs","ro,nosuid,nodev,noexec")
        f.stats[path] = IsolationStat(0,0o644,0,:file)
        @test f.check().uid == 1000
        for (uid,mode,kind) in ((1000,0o644,:file),(0,0o666,:file),(0,0o755,:directory))
            f.stats[path] = IsolationStat(uid,mode,0,kind)
            @test_throws IsolationError f.check()
        end
        f.stats[path] = IsolationStat(0,0o644,0,:file)
        f.files["/proc/self/mountinfo"] = replace(f.files["/proc/self/mountinfo"],
            "$path ro,nosuid,nodev,noexec - tmpfs none ro,nosuid,nodev,noexec"=>
            "$path rw,nosuid,nodev,noexec - tmpfs none rw,nosuid,nodev,noexec")
        @test_throws IsolationError f.check()
    end
end
