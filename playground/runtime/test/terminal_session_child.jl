# Finite, trusted native test child. No container-isolation claim is made.
ccall(:alarm,Cuint,(Cuint,),60)
ccall(:isatty,Cint,(Cint,),0)==1 || exit(70)
ccall(:ioctl,Cint,(Cint,Culong,Cint),0,0x540e,0)==0 || exit(71)
get(ENV,"LCM_TEST_READY","yes")=="yes" &&
    include(joinpath(@__DIR__,"..","..","worker","containers","terminal-ready.jl"))
