using Pkg
Pkg.activate(".")

using Distributed

if ARGS[1] == "--help" || ARGS[1] == "-h" || ARGS[1] == "help"
    println("\n\tFormat: julia tsp_generator.jl dirname n number [threads]

n is the number of cities, number is number of instances. dirname must not
exists. n should be small. Note (n-1)^2 variables are present for QUBO.
threads is the number of threads used, defaults to 1

Generated QUBO and HOBO assume that first cite is visted at time 1.

Following modules are needed: BinaryOptimization, NPZ, LinearAlgebra, Distributed")
    exit()
end

dir_out = ARGS[1]


n = parse(Int, ARGS[2])
@assert n >= 1

repeating = parse(Int, ARGS[3])
@assert repeating >= 1

if isdir(dir_out) || isfile(dir_out)
    error("$dir_out exists, please change or remove it")
    exit(1)
else
    mkdir(dir_out)
end

threads_no = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 0
addprocs(threads_no)
@everywhere using Pkg; 
@everywhere Pkg.activate(".")

@everywhere using LinearAlgebra
@everywhere using BinaryOptimization
@everywhere using NPZ
@everywhere include("tsp_encodings.jl")

@everywhere function simplified_tsp_qubo(w::Matrix{<:Real})
    qubo = tsp_qubo(w, A=2*maximum(w))
    n, _ = indexcapacities(vars(qubo))
    set_one!(qubo, [1,n])
    for i=1:(n-1)
        set_zero!(qubo, [1,i])
        set_zero!(qubo, [i+1,n])
    end
    qubo
end

@everywhere function simplified_tsp_qubo_xy(w::Matrix{<:Real})
    qubo = tsp_qubo(w, Variables("TSP", collect(size(w)), ["timestep", "city"]), 1., 0., 2 * maximum(w))
    n, _ = indexcapacities(vars(qubo))
    set_one!(qubo, [1,n])
    for i=1:(n-1)
        set_zero!(qubo, [1,i])
        set_zero!(qubo, [i+1,n])
    end
    qubo
end

@everywhere function _save_dict_ising(d::Dict, filename::String)
    result = zeros(length(d), 5)

    for (i, (k, v)) = enumerate(collect(d))
        if length(k) >= 1
            result[i, 1] = k[1][1]
            result[i, 2] = k[1][2]
            if length(k) == 2
                result[i, 3] = k[2][1]
                result[i, 4] = k[2][2]
            end
        end
        result[i, 5] = v
    end
    npzwrite(filename, result)
    result
end

@everywhere function save_dict_ising(hobo::HOBO, filename::String)
    save_dict_ising(Ising(hobo), filename)
end

@everywhere function save_dict_ising(ising::Ising, filename::String)
    _save_dict_ising(ising.values, filename)
end

println("arguments are fine. Generating data started ($threads_no threads, $repeating samples)")
@sync @distributed for m = 1:repeating
    println("> $m ($repeating)")
    w = rand(1:10, (n,n))
    qubo_xy = simplified_tsp_qubo_xy(w)  
    save_dict_ising(qubo_xy, "$dir_out/tsp_dict_xy_$m.npz")

    qubo = bom_to_hamiltonian(simplified_tsp_qubo(w))
    npzwrite("$dir_out/qubo_$m.npz", qubo)
end

rmprocs()
