# Fault-injection child: no broker, numerical package or shell subprocess.
mode = only(ARGS)
println("@LCM_EXECUTOR_FRAME@{\"type\":\"bootstrapped\"}")
flush(stdout)
if mode == "malformed"
    readline(stdin)
    println("@LCM_EXECUTOR_FRAME@{\"result\":{}}")
    flush(stdout)
end
while true
    sleep(0.1) # Deliberately never consume the scientific command or shutdown.
end
