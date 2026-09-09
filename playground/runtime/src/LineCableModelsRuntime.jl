"""
    LineCableModelsRuntime

Own application identities, registrations and supervised resource lifetimes.
This package imports neither Bonito nor scientific implementations.
"""
module LineCableModelsRuntime

using Dates, Random, SHA, Sockets, TOML, UUIDs
import FileWatching
using RequiredInterfaces
import AWS, DBInterface, HTTP, JSON3, SQLite, URIs
import LineCableModelsPlaygroundProtocol as Protocol
import LineCableModelsExecutionCore as ExecutionCore

include("Identity.jl")
include("Configuration.jl")
include("Applications.jl")
include("Profiles.jl")
include("OwnedCommands.jl")
include("TerminalProcess.jl")
include(joinpath(@__DIR__, "..", "..", "common", "container_engine.jl"))
include("ContainerHost.jl")
include("EnvironmentFingerprint.jl")
include("DefaultApplications.jl")
include("Store.jl")
include("WorkerStore.jl")
include("LeaseStore.jl")
include("JobStore.jl")
include("JobCancellationStore.jl")
include("WorkerInventory.jl")
include("Assignments.jl")
include("AgentLeases.jl")
include("BrokerPolicy.jl")
include("BrokerTransport.jl")
include("TerminalTransport.jl")
include("LeaseCoordinator.jl")
include("AssignedJobs.jl")
include("RuntimeArtifacts.jl")
include("ControlConfiguration.jl")
include("ControlEvents.jl")
include("ScientificCoordinator.jl")
include("TerminalSockets.jl")
include("TerminalCoordinator.jl")
include("JobCoordinator.jl")
include("ControlService.jl")
include("AgentConfiguration.jl")
include("AgentResources.jl")
include("ResourceJournal.jl")
include("ContainerRecovery.jl")
include("ContainerPolicy.jl")
include("ManagedAgent.jl")
include("NativeRecovery.jl")
include("NativePolicy.jl")
include("ScientificResources.jl")
include("ManagedScientificDriver.jl")
include("ManagedTerminalDriver.jl")
include("TerminalResources.jl")
include("ManagedAgentResources.jl")
include("AgentTerminals.jl")
include("AgentScience.jl")
include("AgentJobs.jl")
include("AgentService.jl")
include("UIHosts.jl")
include("ControlAPI.jl")
include("PublishedSite.jl")
include("Gateway.jl")
include("TerminalGateway.jl")
include("RunSurface.jl")
include("CLI.jl")
include("StartupCompilation.jl")

export AccessDenied, Principal, AbstractIdentityPolicy, ProxyIdentity,
    LocalIdentity, authenticate, authorize_request, RuntimeConfig, read_config,
    RunLimits, ApplicationDefinition, RuntimeRequirement, AbstractApplication,
    HostContext, LocalApplication, ApplicationRegistry, register!, describe,
    ui_command, RuntimeStore, RunRecord, CapacityUnavailable, get_run, list_runs,
    reserve_run!, transition_run!, UIHostSupervisor, start_ui!, stop_ui!, start_gateway,
    ResourceBudget, ProfileDefinition, ProfileRegistry, validate_requirements,
    default_applications, PublishedSite, serve_runtime, runtime_cli

export WorkerTrust, WorkerRegistration, enroll_worker!, list_registrations,
    set_registration_state!, migrate_runtime!

export WorkerInventory, WorkerPresence, probe_worker!, accept_report!, worker_inventory
export CoordinatorIdentity, WorkerIdentity, BrokerEndpoint, BrokerControl,
    broker_permissions, broker_user_config, send_control!, poll_control!
export AutomaticPlacement, PinnedPlacement, DedicatedPlacement, AssignmentLimits,
    LeaseRecord, get_assignment, list_assignments
export AssignmentManager, reserve_assignment!
export JobRecord, get_job, list_jobs, prior_job_submission, reserve_job!, transition_job!
export JobCancellation, job_cancellation, request_job_cancellation!
export AgentLeaseLedger, receive_probe!, handle_lease_control!, expire_agent_leases!,
    complete_agent_cleanup!, agent_lease_usable
export LeaseCoordinator, grant_assignment!, renew_assignment!, release_assignment!,
    accept_lease_ack!, assignment_usable, tick_leases!, reconcile_worker_report!
export BrokerJobs, ensure_worker_streams!, publish_assigned_job!, poll_assigned_job!,
    assigned_result, persist_assigned_result!, AssignedDelivery
export ControlConfig, read_control_config
export LocalRuntimeArtifacts, S3RuntimeArtifacts, ArtifactUnavailable
export ControlService, start_control!, tick_control!, ControlEvents, control_events
export AgentConfig, read_agent_config, AbstractAgentResources, installed_profiles,
    recover_owned!, release_owned!, AgentService, start_agent!, tick_agent!
export EnvironmentFingerprint, native_environment_fingerprint, verify_native_environment
export AbstractScientificDriver, ScientificResources, ScientificOutput, bind_agent!, executor_for!, verify_executor!,
    prepare_assigned!, execute_assigned!, refresh_preparation!, cancel_assigned!, scientific_status, prepared_execution
export cancel_job!
export CommandRunner, CommandResult, CommandFailure, run_owned_command!,
    ContainerEngine, ContainerHostCheck, check_container_host
export ResourceJournal, ResourceReceipt, resource_receipts, resource_name,
    reserve_resource!, bind_resource!, forget_resource!, resource_labels, matches_resource
export container_scope, remove_owned_container!, recover_containers!
export ContainerPolicy, container_create_arguments, create_owned_container!
export ManagedAgentIdentity, agent_service_unit, verify_managed_agent, recover_agent_resources!
export ManagedScientificDriver, ContainerScientificDriver, serve_agent
export ManagedResourceDriver, ManagedScientificView, ManagedTerminalView, AbstractTerminalDriver,
    verify_terminal!, terminal_for!, start_owned_terminal!
export TerminalResources, TerminalSessionLimits, open_terminal!, terminal_status,
    terminal_ready_marker, write_terminal!, read_terminal, resize_terminal!, disconnect_terminal!
export ManagedAgentResources
export BrokerTerminal, watch_terminal!, unwatch_terminal!, send_terminal!, poll_terminal!
export AgentTerminalService, start_terminals!, stop_terminals!, tick_terminals!,
    receive_terminal_command!, stop_terminal!, restart_terminal!
export TerminalCoordinator, request_terminal!, accept_terminal_report!, terminal_flight
export NativeUnitIdentity, native_scope, native_resource_description,
    inspect_native_unit, remove_owned_native!, recover_native!
export NativeHostCheck, check_native_host, NativePolicy, native_launch_command, verify_native_service!
export AgentScientificService, receive_scientific!, start_science!, tick_science!, stop_science!
export ScientificCoordinator, request_scientific!, accept_scientific_report!, remote_scientific_status
export AgentJobService, start_jobs!, tick_jobs!, stop_jobs!, accept_assigned_delivery!, progress_assigned_delivery!
export JobCoordinator, submit_job!, owned_job_result

end
