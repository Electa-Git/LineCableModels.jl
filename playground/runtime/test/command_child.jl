mode = only(ARGS)
if mode == "echo"
    print(stdout, "private-stdout")
    print(stderr, "private-stderr")
    exit(17)
elseif mode == "stdin"
    print(eof(stdin) ? "closed" : "open")
elseif mode == "flood"
    while true
        print(stdout, repeat("private-flood", 1000))
        flush(stdout)
        print(stderr, repeat("private-flood", 1000))
        flush(stderr)
    end
elseif mode == "wait"
    sleep(30)
else
    exit(71)
end
