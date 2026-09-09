# Loaded by the successful container guard, not a replacement for that guard.
# Keep this fixed hook separate so the real Julia REPL startup path can also be
# tested without pretending a native test process is a sandboxed container.
let token=get(ENV,"LCM_TERMINAL_READY","")
    occursin(r"^[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}$",token) || exit(78)
    atreplinit() do repl
        write(stdout,"\x1elcm-terminal-ready:",token,"\x1f")
        flush(stdout)
    end
end
