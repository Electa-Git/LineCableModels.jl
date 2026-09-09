using Test, UUIDs, LineCableModelsExecutionCore

function native_isolation_fixture()
    container = isolation_fixture()
    group = "/user.slice/user-1000.slice/user@1000.service/app.slice/lcm-exec-$(uuid4()).service"
    identity = NativeIdentity(1000,1000,group)
    files = copy(container.files)
    files["/proc/self/cgroup"] = "0::" * group * "\n"
    for name in ("cpu.max","memory.max","memory.swap.max","pids.max")
        files["/sys/fs/cgroup" * group * "/" * name] = pop!(files,"/sys/fs/cgroup/" * name)
    end
    mount = container.mount
    files["/proc/self/mountinfo"] = join([
        mount("/","xfs","ro"),mount("/home","ext4","ro"),mount("/proc","proc","rw,nosuid,nodev,noexec"),
        mount("/sys","sysfs","ro"),mount("/sys/fs/cgroup","cgroup2","ro"),mount("/dev","tmpfs","ro"),
        mount("/dev/pts","devpts","rw,nosuid,noexec"),mount("/tmp","tmpfs","rw,nosuid,nodev,noexec"),
        mount("/dev/shm","tmpfs","rw,nosuid,nodev,noexec")],'\n')
    disks = Dict{String,Any}(container.disks)
    disks["/sys/fs/cgroup"] = (ftype=0x63677270,)
    check() = E.check_native_isolation(container.limits,identity,p->files[p],p->disks[p])
    return (;files,disks,identity,check,mount,limits=container.limits,group)
end

@testset "native entry shares quotas but requires exact native identity" begin
    @test ExecutorLimits === ContainerLimits
    f = native_isolation_fixture()
    @test f.check() == (schema_version=1,cpus=1.0,memory_bytes=1024^3,pids=128,scratch_bytes=256*1024^2,uid=1000)
    @test !occursin(f.group,repr(MIME"text/plain"(),f.identity))
    for (uid,gid,group) in ((0,1000,f.group),(true,1000,f.group),(1000,0,f.group),
            (1001,1000,f.group),(1000,1000,f.group*"/child"),(1000,1000,"/host"),
            (1000,1000,replace(f.group,r"lcm-exec-.*"=>"lcm-exec-"*repeat("-",36)*".service")))
        @test_throws ArgumentError NativeIdentity(uid,gid,group)
    end
    for (name,value) in (("cpu.max","max 100000"),("cpu.max","100001 100000"),("cpu.max","1 1"),
            ("memory.max","max"),("memory.max",string(1024^3+1)),("memory.swap.max","1"),
            ("pids.max","max"),("pids.max","129"))
        f = native_isolation_fixture(); f.files["/sys/fs/cgroup"*f.group*"/"*name] = value
        @test_throws IsolationError f.check()
    end
    for (old,new) in (("Uid:\t1000 1000 1000 1000","Uid:\t0 0 0 0"),
            ("Gid:\t1000 1000 1000 1000","Gid:\t0 0 0 0"),
            ("CapBnd:\t0000000000000000","CapBnd:\t0000000000000001"),("NoNewPrivs:\t1","NoNewPrivs:\t0"))
        f = native_isolation_fixture(); f.files["/proc/self/status"] = replace(f.files["/proc/self/status"],old=>new)
        @test_throws IsolationError f.check()
    end
    f = native_isolation_fixture(); f.files["/proc/self/cgroup"] *= "extra"
    @test_throws IsolationError f.check()
    f = native_isolation_fixture(); f.disks["/sys/fs/cgroup"] = (ftype=0x1234,)
    @test_throws IsolationError f.check()
    f = native_isolation_fixture(); f.files["/proc/net/dev"] *= "eth0: 0\n"
    @test_throws IsolationError f.check()
    f = native_isolation_fixture(); f.files["/proc/self/limits"] = replace(f.files["/proc/self/limits"],"1024 1024"=>"1024 unlimited")
    @test_throws IsolationError f.check()
    for path in ("/","/home","/sys","/sys/fs/cgroup","/dev")
        f = native_isolation_fixture()
        f.files["/proc/self/mountinfo"] = replace(f.files["/proc/self/mountinfo"],Regex("( / " * path * " )ro( - [^\\n]+ )ro")=>s"\1rw\2rw")
        @test_throws IsolationError f.check()
    end
    for (path,fs,options) in (("/data","ext4","rw"),("/tmp/nested","tmpfs","rw"),
            ("/proc/private","tmpfs","rw"),("/tmp","tmpfs","rw,nosuid,nodev,noexec"))
        f = native_isolation_fixture(); f.files["/proc/self/mountinfo"] *= "\n" * f.mount(path,fs,options)
        @test_throws IsolationError f.check()
    end
    f = native_isolation_fixture(); f.disks["/dev/shm"] = (ftype=0x01021994,bsize=4096,blocks=17)
    @test_throws IsolationError f.check()
    f = native_isolation_fixture(); f.disks["/tmp"] = (ftype=0x01021994,bsize=4096,blocks=65536)
    @test_throws IsolationError f.check()
    f = native_isolation_fixture(); f.files["/proc/self/mountinfo"] = replace(f.files["/proc/self/mountinfo"],"rw,nosuid,nodev,noexec"=>"rw,nosuid,nodev")
    @test_throws IsolationError f.check()
end
