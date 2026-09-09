@testset "terminal frames are private, bounded and assignment-fenced" begin
    fence=AssignmentFence(string(uuid4()),string(uuid4()),"alice","terminal","worker-a",
        string(uuid4()),string(uuid4()),"julia-terminal","1.0.0",repeat("a",64),1)
    writer,session=string(uuid4()),string(uuid4())
    open=TerminalCommand("2.0",string(uuid4()),fence,1,"open",nothing,writer,0,0,100,30,UInt8[])
    input=replace_record(open;action="input",session_id=session,input_sequence=1,columns=0,rows=0,
        bytes=collect(codeunits("private-code-λ\r")))
    status=replace_record(input;action="status",writer_id=nothing,input_sequence=0,bytes=UInt8[])
    read=replace_record(status;action="read",after=42)
    report=TerminalReport("2.0",read.request_id,fence,1,true,"accepted",session,"ready",true,1,45,48,false,
        UInt8[0x00,0x80,0xff],nothing,false)
    for record in (open,input,status,read,report)
        @test decode_terminal_message(typeof(record),encode_message(record))==record
        @test !occursin(writer,repr(record)) && !occursin("private-code",repr(record))
    end
    for action in ("resize","keepalive","disconnect","stop","restart")
        command=replace_record(input;action,input_sequence=0,bytes=UInt8[],
            columns=action in ("resize","restart") ? 80 : 0,rows=action in ("resize","restart") ? 24 : 0)
        @test decode_terminal_message(TerminalCommand,encode_message(command))==command
    end
    for changes in ((action="eval",),(session_id=session,),(writer_id=nothing,),(revision=0,),
            (columns=1001,),(rows=0,),(after=1,),(input_sequence=1,),(bytes=UInt8[1],),(protocol_version="1.0",))
        @test_throws ArgumentError validate(replace_record(open;changes...))
    end
    @test_throws ArgumentError validate(replace_record(input;session_id=nothing))
    @test_throws ArgumentError validate(replace_record(input;input_sequence=0))
    @test_throws ArgumentError validate(replace_record(input;bytes=zeros(UInt8,8193)))
    @test_throws ArgumentError validate(replace_record(status;writer_id=writer))
    @test_throws ArgumentError validate(replace_record(status;action="keepalive"))
    @test_throws ArgumentError validate(replace_record(read;after=-1))
    for changes in ((accepted=false,),(phase="installed",),(reason="/private/path",),
            (cursor=2,),(output_sequence=1,),(failure="private exception",),(session_id=nothing,),
            (bytes=zeros(UInt8,8193),))
        @test_throws ArgumentError validate(replace_record(report;changes...))
    end
    for (field,value) in (("revision",true),("input_sequence",true),("columns",1.5),
            ("bytes",[256]),("bytes",[-1]),("bytes",[true]),("bytes","code"),("command","sh"))
        data=JSON3.read(encode_message(input),Dict{String,Any});data[field]=value
        @test_throws ArgumentError decode_terminal_message(TerminalCommand,JSON3.write(data))
    end
    @test_throws ArgumentError decode_terminal_message(TerminalCommand,fill(UInt8(' '),MAX_TERMINAL_FRAME_BYTES+1))
    subject=terminal_subject(fence,:command)
    @test startswith(subject,"lcm.terminal.v2.worker-a.") && !startswith(subject,"lcm.jobs.")
    for changed in (replace_record(fence;worker_id="worker-b"),replace_record(fence;worker_boot=string(uuid4())),
            replace_record(fence;lease_id=string(uuid4())),replace_record(fence;generation=2))
        @test terminal_subject(changed,:command)!=subject
    end
    @test terminal_subject(fence,:report)!=subject
    @test_throws ArgumentError terminal_subject(fence,:wildcard)
end
