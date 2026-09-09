### tls.jl
#
# Copyright (C) 2023 Jakub Wronowski.
#
# Maintainer: Jakub Wronowski <jakubwro@users.noreply.github.com>
# Keywords: nats, nats-client, julia
#
# This file is a part of NATS.jl.
#
# License is MIT.
#
### Commentary:
#
# This file contains utilities for handling TLS handshake.
#
### Code:

function upgrade_to_tls(
        sock::Sockets.TCPSocket,
        server_name::AbstractString,
        ca_cert_path::Union{String, Nothing},
        client_cert_path::Union{String, Nothing},
        client_key_path::Union{String, Nothing},
)
    entropy = MbedTLS.Entropy()
    rng = MbedTLS.CtrDrbg()
    MbedTLS.seed!(rng, entropy)
    ctx = MbedTLS.SSLContext()
    conf = MbedTLS.SSLConfig()
    MbedTLS.config_defaults!(conf)
    MbedTLS.authmode!(conf, MbedTLS.MBEDTLS_SSL_VERIFY_REQUIRED)
    MbedTLS.rng!(conf, rng)

    # function show_debug(level, filename, number, msg)
    #     @show level, filename, number, msg
    # end
    
    # MbedTLS.dbg!(conf, show_debug)
    
    if !isnothing(ca_cert_path)
        MbedTLS.ca_chain!(conf, MbedTLS.crt_parse_file(ca_cert_path))
    end

    if !isnothing(client_cert_path) && !isnothing(client_key_path)
        cert = MbedTLS.crt_parse_file(client_cert_path)
        key = MbedTLS.parse_keyfile(client_key_path)
        MbedTLS.own_cert!(conf, cert, key)
    elseif !isnothing(client_cert_path) || !isnothing(client_key_path)
        error("Both the TLS client certificate and key must be provided.")
    end

    MbedTLS.setup!(ctx, conf)
    MbedTLS.hostname!(ctx, server_name)
    MbedTLS.set_bio!(ctx, sock)

    MbedTLS.handshake(ctx)

    get_tls_input_buffered(ctx), ctx
end

function copy_tls_input!(io, ssl)
    try
        while !eof(ssl)
            write(io, readavailable(ssl))
        end
    catch error
        # MbedTLS.eof may resume after another task closed its read side.
        # Match only that dependency's closed-reader guard, not TLS alerts,
        # verification failures or arbitrary errors on a closed connection.
        error isa Base.IOError && error.code == 0 &&
            error.msg == "`ssl_unsafe_read` requires `isreadable(::SSLContext)`" &&
            !isreadable(ssl) || rethrow()
    finally
        close(io)
    end
    return nothing
end

function get_tls_input_buffered(ssl)
    io = Base.BufferStream()
    t = Threads.@spawn :interactive disable_sigint() do
        copy_tls_input!(io, ssl)
    end
    errormonitor(t)
    BufferedInputStream(io, 1)
end
