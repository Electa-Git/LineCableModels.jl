using LineCableModelsRuntime, UUIDs
const R = LineCableModelsRuntime
journal = ResourceJournal(only(ARGS), "worker-a")
try
    fence = R.Protocol.AssignmentFence(string(uuid4()), string(uuid4()), "alice", "main", "worker-a",
        string(uuid4()), string(uuid4()), "line-parameters", "1.0.0", repeat("a",64), 1)
    reserve_resource!(journal, fence, :podman, repeat("e",64))
    println("journal-ready")
    flush(stdout)
    readline(stdin)
finally
    close(journal)
end
