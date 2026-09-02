struct WorkerStatus
    client::BrokerClient
    operation::String
end

function Bonito.jsrender(session::Session, status::WorkerStatus)
    label = map(
        session,
        status.client.connection_status,
        status.client.available_capabilities
    ) do connection, workers
        connection == :online || return "Broker: $(connection)"
        compatible = count(
            capabilities -> status.operation in capabilities,
            values(workers)
        )
        compatible == 0 && return "No compatible worker"
        return "$compatible compatible worker$(compatible == 1 ? "" : "s")"
    end
    state_class = map(session, status.client.connection_status) do state
        "lc-worker-state lc-worker-state-$state"
    end
    return Bonito.jsrender(
        session,
        DOM.span(label; class=state_class, role="status")
    )
end
