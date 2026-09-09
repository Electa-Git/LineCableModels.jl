# Finite trusted fixture for PTY behavior, never a production terminal backend.
# An independent kernel deadline bounds even a broken parent test supervisor.
ccall(:alarm,Cuint,(Cuint,),60)
using REPL
ccall(:isatty,Cint,(Cint,),0)==1 && ccall(:isatty,Cint,(Cint,),1)==1 || exit(70)
ccall(:ioctl,Cint,(Cint,Culong,Cint),0,0x540e,0)==0 || exit(71) # Linux TIOCSCTTY.
terminal=REPL.Terminals.TTYTerminal("xterm-256color",stdin,stdout,stderr)
mode=only(ARGS)
if mode=="repl"
    Base.exit_on_sigint(false) # Match Julia --interactive, rather than script SIGINT exit.
    println("REPL-FIXTURE");flush(stdout)
    REPL.run_repl(REPL.LineEditREPL(terminal,true))
else
    REPL.Terminals.raw!(terminal,true)
    println("READY");flush(stdout)
    if mode=="echo"
        while !eof(stdin)
            line=readline(stdin)
            line=="quit" && break
            println("ECHO:",line);flush(stdout)
        end
    elseif mode=="stall"
        sleep(30)
    elseif mode=="flood"
        block=repeat("private-terminal-output-",1024)
        for _ in 1:256
            write(stdout,block);flush(stdout)
        end
    else
        exit(72)
    end
end
