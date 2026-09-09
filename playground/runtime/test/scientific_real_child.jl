ccall(:getppid, Cint, ()) == parse(Int, ENV["LCM_EXECUTOR_PARENT_PID"]) || exit(70)
# Only the explicit test driver chooses these two installed child environments.
if only(ARGS) == "line-parameters"
    using LineCableModelsLineParameters
    LineCableModelsLineParameters.main()
elseif only(ARGS) == "power-flow"
    using LineCableModelsPowerFlow
    LineCableModelsPowerFlow.main()
else
    error("unregistered scientific test profile")
end
