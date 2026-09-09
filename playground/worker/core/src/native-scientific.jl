# Load only the approved package identity after the fixed kernel guard passed.
# Package UUID/name come from read-only source fingerprint verification.
isdefined(Main,:LCM_NATIVE_ISOLATION) || exit(78)
length(ARGS) == 2 || exit(78)
const LCM_NATIVE_PACKAGE = Base.require(Base.PkgId(Base.UUID(ARGS[1]),ARGS[2]))
empty!(ARGS)
Base.invokelatest(getfield(LCM_NATIVE_PACKAGE,:main))
