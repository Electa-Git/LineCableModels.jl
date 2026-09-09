# Deterministically reproduce the dependency's eof/close interleaving without
# weakening certificate verification or needing a probabilistic network race.
struct TLSReaderFailure <: IO
    readable::Bool
    error::Exception
end
Base.isreadable(io::TLSReaderFailure) = io.readable
Base.eof(io::TLSReaderFailure) = throw(io.error)

@testset "TLS reader closes cleanly without swallowing transport failures" begin
    closed_guard = Base.IOError("`ssl_unsafe_read` requires `isreadable(::SSLContext)`", 0)
    io = Base.BufferStream()
    @test RT.NATS.copy_tls_input!(io, TLSReaderFailure(false, closed_guard)) === nothing
    @test !isopen(io) && eof(io)
    for source in (TLSReaderFailure(true, closed_guard),
            TLSReaderFailure(false, Base.IOError(closed_guard.msg, 1)),
            TLSReaderFailure(false, Base.IOError("TLS alert", 0)),
            TLSReaderFailure(false, ErrorException("unexpected reader failure")))
        output = Base.BufferStream()
        caught = try
            RT.NATS.copy_tls_input!(output, source)
        catch error
            error
        end
        @test caught === source.error
        @test !isopen(output) && eof(output)
    end
    output = Base.BufferStream()
    @test RT.NATS.copy_tls_input!(output, IOBuffer("bounded data")) === nothing
    @test !isopen(output)
    @test String(read(output)) == "bounded data"
end
