# Keep first-use JSON compilation outside live acknowledgement deadlines.
# These declarations compile passive codecs; they send no messages, create no
# authority and neither connect to a broker nor load a numerical environment.
for record in (WorkerProbe, WorkerAnnouncement, LeaseControl, LeaseAcknowledgement,
        ScientificCommand, ScientificReport, AssignedJob, AssignedResult)
    precompile(encode_message, (record,))
    precompile(decode_runtime_message, (Type{record}, String))
    precompile(decode_runtime_message, (Type{record}, Vector{UInt8}))
end
