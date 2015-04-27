using BinDeps
using Compat

vers = "20U1"

tagfile = "installed_vers"
target = "libbid$WORD_SIZE.$(Libdl.dlext)"
if !isfile(tagfile) || !isfile(target) || readchomp(tagfile) != "$vers $WORD_SIZE"
    if OS_NAME == :Windows
        error("no binary DLL available yet for Windows")
#        run(download_cmd("http://ab-initio.mit.edu/dfp/libbid$WORD_SIZE-$vers.dll", target))
#    elseif OS_NAME == :Darwin
#        run(download_cmd("http://ab-initio.mit.edu/dfp/libbid$WORD_SIZE-$vers.dylib", target))
    else
        tarball = "IntelRDFPMathLib$vers.tar.gz"
        srcdir = "IntelRDFPMathLib$vers/LIBRARY"
        if !isfile(tarball)
            run(download_cmd("http://www.netlib.org/misc/intel/$tarball", tarball))
        end
        run(unpack_cmd(tarball, ".", ".gz", ".tar"))
        cd(srcdir) do
            println("COMPILING LIBBID...")
            open("makefile.shared", "w") do f
                println(f, "include makefile\n\n../../$target: \$(ALL_BID_OBJS)\n\t\$(CC) -shared -o \$@ \$(ALL_BID_OBJS)\n")
            end
            run(`make -f makefile.shared CC=gcc CFLAGS_OPT=-O2 CFLAGS_CC=-fPIC CALL_BY_REF=0 GLOBAL_RND=1 GLOBAL_FLAGS=1 UNCHANGED_BINARY_FLAGS=1 ../../$target`)
        end
    end
    open(tagfile, "w") do f
        println(f, "$vers $WORD_SIZE")
    end
end
