"""
    LineCableModelsExecutionCore

Share operation validation, scientific-process framing, cancellation and caches
between the legacy worker and profile-specific executors. This package imports
neither numerical packages nor broker, artifact-client, or UI implementations.
"""
module LineCableModelsExecutionCore

using Dates, JSON3, Logging, SHA, UUIDs
using LineCableModelsPlaygroundProtocol
using RequiredInterfaces

# These are the existing implementation sources, not parallel copies. Keep the
# monorepo worker tree together when installing a profile or building an image.
include("../../src/OperationRegistry.jl")
include("../../src/Cache.jl")
include("../../src/Executor.jl")
include("Profiles.jl")
include("ContainerIsolation.jl")
include("NativeIsolation.jl")

export OperationSpec, OperationRegistry, register!, registered_operation,
    PermanentOperationError, RetryableOperationError, required, bounded_real,
    bounded_integer, input_hash, normalize_wire, ExecutionContext, CancellationToken,
    cancel!, check_canceled, progress!, log!, warning!, execute_operation,
    PreparedResourceCache, prepare_resource!, prepared_status, clear_prepared!,
    ExecutorSupervisor, start_executor!, stop_executor!, execute_supervised!, prepare_supervised!, inspect_preparation!,
    AbstractScientificProfile, PreparedWorkload, operation_registry, load_profile!, validate_preparation, prepare!,
    cleanup!, profile_executor_main
export ContainerLimits, IsolationError, verify_container_isolation
export ExecutorLimits, NativeIdentity, verify_native_isolation

end
