@testsnippet TemporaryFEMRuntime begin
    # The runner restores each item's working directory; the item owns cleanup.
    const fem_test_working_directory=pwd()
    const fem_test_runtime_directory=mktempdir()
end
