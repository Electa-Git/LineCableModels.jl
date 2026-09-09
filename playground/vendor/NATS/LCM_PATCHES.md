# LineCableModels NATS transport patch

This directory is an auditable local fork of the MIT-licensed registered
`NATS.jl` 0.1.1 source (`git-tree-sha1`
`a1cdf34ba90ee5cd2658e487d3277ffafee712ce`). The upstream repository is
<https://github.com/jakubwro/NATS.jl>.

The fork changes connection and request/reply transport:

1. accept a mandatory-TLS server advertising `tls_required` without requiring
   the obsolete optional `tls_available` flag;
2. configure the client certificate before constructing the TLS context and
   reject incomplete certificate/key pairs;
3. set and verify the expected server hostname through MbedTLS;
4. support an explicit `NATS_TLS_SERVER_NAME` override for connections whose
   network address differs from the certificate identity;
5. treat an intentional drain during initial reconnect as normal cleanup,
   avoiding false error logs when an offline publisher shuts down;
6. bound each TCP/INFO/TLS/authentication handshake with `connect_timeout`
   (default 10 seconds). An owned timer closes only the attempt's socket;
   success/failure cancels it. Initial connection and reconnect share this path;
7. honor the client's TLS requirement when a server advertises optional TLS.
   An optional server must not cause plaintext CONNECT credentials when the
   client explicitly requires encryption;
8. honor `send_enqueue_when_disconnected=false` during CONNECTING as well as
   DISCONNECTED. Lease renewals and terminal bytes must not silently queue for
   later replay while a connection is being re-established;
9. support a validated per-connection `inbox_prefix` (default `inbox.`), including
   JetStream requests, and size request reply buffers to the requested reply
   count. Private worker credentials no longer need a shared wildcard inbox;
10. forward the requested JetStream acknowledgement body. Previously `-NAK`,
    `+TERM` and `+WPI` were silently sent as an empty positive acknowledgement;
11. omit raw connection URLs from reconnect diagnostics and malformed-host
    errors. Legacy URLs can contain passwords in their user-info component;
12. clean subscription resources when SUB/UNSUB transmission fails. A bounded
    reconnect policy must not leave reply channels behind after failed requests;
13. discard the uncertain send buffer on reconnect when disconnected enqueue is
    disabled. Subscriptions are reconstructed, but already queued terminal or
    control bytes must not be replayed onto a replacement connection;
14. close request timers and reply subscriptions in a finally block, including
    publication failure and early error responses. Reply expiry closes its queue
    immediately instead of waiting for an application reader to drain it;
15. end the buffered TLS reader normally for MbedTLS's specific closed-reader
    guard when its read side actually closed during `eof`. Other IO errors,
    TLS alerts, verification failures and errors on a still-readable context
    remain errors. The consumer buffer closes on every exit path. The TLS copy
    loop and buffered receiver have passive precompile declarations so their
    first specialization does not consume the initial PONG deadline.

16. compile the configured connection controller before returning the connection.
    A measured 2.282-second first specialization otherwise delayed sender/receiver
    startup beyond the caller's 2-second PONG deadline. The compile request sends
    nothing and changes no handshake/lease timeout; failure closes the owned
    socket and finishes the drain waiter.

The runtime integration suite tests peers that accept TCP but never send INFO,
and peers that send INFO but never complete authentication. TLS verification
still follows the same verified transport path; this adds no insecure fallback.

The application pins this directory through Julia's `[sources]` mechanism.
No certificate verification bypass is provided.
