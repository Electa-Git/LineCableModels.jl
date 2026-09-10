using TestItemRunner, LineCableModels
TestItemRunner.run_tests(pkgdir(LineCableModels);
    filter = item -> :gauntlet in item.tags || :gauntlet_toolkit in item.tags,
    verbose = true)
