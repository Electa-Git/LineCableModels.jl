temporary() = joinpath(tempdir(), "x")
staged() = mktempdir()
named() = Base.tempname()
home() = homedir()
