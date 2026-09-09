ccall(:getppid, Cint, ()) == parse(Int, ENV["LCM_UI_PARENT_PID"]) || exit(70)
# Never writes readiness; the supervisor must time out and reclaim this run.
sleep(600)
