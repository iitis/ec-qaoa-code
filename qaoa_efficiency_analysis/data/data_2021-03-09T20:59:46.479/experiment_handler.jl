
using SHA
using Dates

function is_important(filename)
    res = false
    res = res || (length(filename)>3 && filename[end-2:end] == ".jl")
    res = res || (length(filename)>5 && filename[end-4:end] == ".toml")
    res = res && !occursin("plot", filename)
    res
end

function check_if_same_files(dir_out::String)
    for filename = readdir()
        if isfile(filename) && is_important(filename)
            println("> $filename under consideration")
            @assert isfile("$dir_out/$filename") "$dir_out/$filename doesn't exist, and it should"
            sha1 = ""
            sha2 = ""
            open("$filename") do f
                sha1 = sha2_256(f)
            end
            open("$dir_out/$filename") do f
                sha2 = sha2_256(f)
            end
            @assert sha1 == sha2 "$filename and $dir_out/$filename doesn't match, and it should - versions are likely different"  
        else
            println("> $filename skipped (not important)")
        end
    end
    nothing
end

function copy_files(dir_out::String)
    println("The directory exists - File checking")
    for filename = readdir()
        if is_important(filename)
            cp(filename, "$dir_out/$filename")
        end
    end
    nothing
end

function save_readme(dir_out)
    open("$dir_out/README_$(Dates.now()).md", "w") do f
        write(f, "$PROGRAM_FILE $ARGS")
    end
    nothing
end

threads_no = parse(Int, ARGS[findfirst(x->x=="-p", ARGS)+1])
dir_in = ARGS[findfirst(x->x=="-in",ARGS)+1]

if findfirst(x->x=="-outold", ARGS) != nothing
    dir_out = ARGS[findfirst(x->x=="-outold", ARGS)+1]
    @assert isdir(dir_out) "$dir_out doesn't exist"

    dir_out *= "/"
    println("Version checking")
    check_if_same_files(dir_out)
    dir_in = "$dir_out/$dir_in"
    dir_out *= "/data"
else
    dir_out = ARGS[findfirst(x->x=="-outnew", ARGS)+1]
    @assert isdir(dir_out) "$dir_out doesn't exist"

    dir_out *= "/data_$(Dates.now())"
    mkdir(dir_out) 
    copy_files(dir_out)
    cp(dir_in, dir_out*"/"*dir_in)
    dir_in = "$dir_out/$dir_in"
    dir_out *= "/data"
    mkdir(dir_out)
end



save_readme(dir_out)
@show dir_out
@show dir_in
@show threads_no

println()