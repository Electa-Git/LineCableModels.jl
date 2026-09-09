using Sockets

# Use the same browser renderer and protected gateway with a real leased REPL.
# Isolation belongs to the caller's driver/host gate, not to this browser check.
function terminal_live_browser(admit, service, supervisor, principal, run, directory)
    site_dir=mkpath(joinpath(directory,"browser-site"))
    config=JSON3.write((kind="terminal",run_id=string(run.id),role="terminal",title="Julia REPL",rows=18))
    css=("brand.css","control-contract.css","forms.css","runtime-controls.css","runtime-terminal.css","runtime-terminal.bundle.css")
    scripts=("runtime-client.js","runtime-terminal-client.js","runtime-terminal.js","runtime-terminal.bundle.js")
    html="""<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
      <title>Live private terminal fixture</title>$(read(RT.RUNTIME_THEME_INIT,String))
      $(join(("<link rel=\"stylesheet\" href=\"/runtime/assets/$name\">" for name in css),"\n"))
      <style>body{margin:12px;background:var(--lc-bg);color:var(--lc-text)}main{max-width:1000px;margin:auto}</style>
      </head><body><main><div id="terminal" data-lcm-runtime-terminal='$(RT.html_text(config))'></div></main>
      <script id="terminal-config" type="application/json">$config</script>
      $(join(("<script src=\"/runtime/assets/$name\"></script>" for name in scripts),"\n"))
      <script>
      const OriginalTerminal=LineCableModelsTerminalVendor.Terminal;
      globalThis.LineCableModelsTerminalVendor={...LineCableModelsTerminalVendor,Terminal:class extends OriginalTerminal {
        constructor(...args){super(...args);globalThis.__liveTerminal=this;}
      }};
      LineCableModelsTerminal.mount(document.querySelector('#terminal'));
      </script></body></html>"""
    write(joinpath(site_dir,"index.html"),html)
    listener=Sockets.listen(ip"127.0.0.1",0)
    port=Int(last(Sockets.getsockname(listener)));close(listener)
    origin="http://127.0.0.1:$port"
    server=start_gateway(supervisor,LocalIdentity(origin,principal);port,site=PublishedSite(site_dir),control=service)
    script=joinpath(@__DIR__,"..","..","test","integration","runtime_terminal_browser.mjs")
    process=nothing
    try
        # This second, local-identity gateway shares the fixture's compiler with
        # both lease peers. Complete its cold HTTP setup before granting browser
        # authority; the first HTTP listener specialization alone takes ~3 s,
        # longer than the unchanged 2 s lease-ack deadline. No terminal is opened
        # by these read-only startup checks.
        @test all(lease->lease.fence.role!="terminal",
            list_assignments(supervisor.store,principal;run_id=run.id))
        for path in ("/health","/runtime/api/control","/runtime/api/runs/$(run.id)/assignments")
            @test HTTP.get(origin*path;proxy=nothing,cookies=false,request_timeout=10).status==200
        end
        admit()
        # Container startup is not a 7.5-second DOM transition. Reuse the
        # shipped lifecycle bound for initial startup and explicit restart;
        # ordinary interactions and reconnect retain their short browser bound.
        startup_seconds=RT.TerminalSessionLimits().startup_seconds
        command=addenv(`node $script`,"LCM_TERMINAL_TEST_ORIGIN"=>origin*"/",
            "LCM_TERMINAL_TEST_STARTUP_MS"=>string(round(Int,1000startup_seconds)))
        process=run_process=Base.run(pipeline(command;stdout=stdout,stderr=stderr);wait=false)
        @test timedwait(()->process_exited(run_process),2startup_seconds+100;pollint=0.1)==:ok
        if !process_exited(process)
            kill(process,Base.SIGTERM)
            timedwait(()->process_exited(process),10)==:ok || kill(process,Base.SIGKILL)
        end
        wait(process)
        if !success(process)
            diagnostics=(control=RT.control_snapshot(service,principal),
                assignments=[RT.assignment_payload(service,principal,lease)
                    for lease in list_assignments(supervisor.store,principal;run_id=run.id)],
                events=RT.control_events(service.events,principal))
            path=joinpath(directory,"terminal-browser-runtime.json")
            write(path,JSON3.write(diagnostics))
            println(stderr,"Terminal runtime diagnostics: ",path)
        end
        @test success(process)
    finally
        if process !== nothing && !process_exited(process)
            kill(process,Base.SIGKILL);wait(process)
        end
        close(server)
    end
end
