# LineCableModels NATS transport patch

This directory is an auditable local fork of the MIT-licensed registered
`NATS.jl` 0.1.1 source (`git-tree-sha1`
`a1cdf34ba90ee5cd2658e487d3277ffafee712ce`). The upstream repository is
<https://github.com/jakubwro/NATS.jl>.

The fork changes only the connection transport:

1. accept a mandatory-TLS server advertising `tls_required` without requiring
   the obsolete optional `tls_available` flag;
2. configure the client certificate before constructing the TLS context and
   reject incomplete certificate/key pairs;
3. set and verify the expected server hostname through MbedTLS;
4. support an explicit `NATS_TLS_SERVER_NAME` override for connections whose
   network address differs from the certificate identity;
5. treat an intentional drain during initial reconnect as normal cleanup,
   avoiding false error logs when an offline publisher shuts down.

The application pins this directory through Julia's `[sources]` mechanism.
No certificate verification bypass is provided.
