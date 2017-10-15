using BinDeps

vers = "20U1"

url="https://bintray.com/artifact/download/julialang/generic/"
tagfile = "installed_vers"
target = "libbid$(Sys.WORD_SIZE).$(Libdl.dlext)"
if !isfile(tagfile) || !isfile(target) || readchomp(tagfile) != "$vers:$(Sys.WORD_SIZE)"
    if Sys.KERNEL == :NT
        # binary for Windows was cross-compiled with mingw using:
        # 32-bit: CC_NAME_INDEX=3 CC_INDEX=3 _HOST_OS=Windows_NT _HOST_ARCH=x86 _NUM_CPUS=1 CC=i686-w64-mingw32-gcc CFLAGS_OPT="-O2 -DBID_THREAD= -DBID_MS_FLAGS" CALL_BY_REF=0 GLOBAL_RND=1 GLOBAL_FLAGS=1 UNCHANGED_BINARY_FLAGS=1
        # 64-bit: CC_NAME_INDEX=3 CC_INDEX=3 _HOST_OS=Windows_NT _HOST_ARCH=x86_64 _NUM_CPUS=1 CC=x86_64-w64-mingw32-gcc CFLAGS_OPT="-O2 -DBID_THREAD= -DBID_MS_FLAGS" CALL_BY_REF=0 GLOBAL_RND=1 GLOBAL_FLAGS=1 UNCHANGED_BINARY_FLAGS=1
        # with the makefile.shared file:
        #   include makefile
        #   libbid.dll: $(ALL_BID_OBJS)
        #           $(CC) -shared -o $@ $(ALL_BID_OBJS)
        run(download_cmd(url*"libbid$(Sys.WORD_SIZE)-$vers.dll", target))
    elseif Sys.KERNEL == :Darwin
        run(download_cmd(url*"libbid$(Sys.WORD_SIZE)-$vers.dylib", target))
    else
        tarball = "IntelRDFPMathLib$vers.tar.gz"
        srcdir = "IntelRDFPMathLib$vers/LIBRARY"
        if !isfile(tarball)
            run(download_cmd(url*"$tarball", tarball))
        end
        run(unpack_cmd(tarball, ".", ".gz", ".tar"))
        cd(srcdir) do
            println("COMPILING LIBBID...")
            open("makefile.shared", "w") do f
                println(f, "include makefile\n\n../../$target: \$(ALL_BID_OBJS)\n\t\$(CC) -shared -o \$@ \$(ALL_BID_OBJS)\n")
            end
            run(`make -f makefile.shared CC=gcc CFLAGS_OPT="-O2 -fPIC" CALL_BY_REF=0 GLOBAL_RND=1 GLOBAL_FLAGS=1 UNCHANGED_BINARY_FLAGS=1 ../../$target`)
        end
    end
    open(tagfile, "w") do f
        println(f, "$vers:$(Sys.WORD_SIZE)")
    end
end
