using Pkg
Pkg.activate(".")
using Distributed

if ARGS[1] == "--help" || ARGS[1] == "-h" || ARGS[1] == "help"
    println("\n\tFormat: julia hamilton_generator.jl dirname nmax [threads]

nmax is the maximal number of cities (procedure starts always from three).
dirname must not exists. nmax should be small. Note (n-1)^2 variables are present for QUBO.
threads is the number of threads used, defaults to 1

Generated QUBO and HOBO assume that first cite is visted at time 1.

Following modules are needed: BinaryOptimization, NPZ, LinearAlgebra, Distributed")
    exit()
end

dir_out = ARGS[1]


nmax = parse(Int, ARGS[2])
@assert nmax >= 3

if isdir(dir_out) || isfile(dir_out)
    error("$dir_out exists, please change or remove it")
    exit(1)
else
    mkdir(dir_out)
end

threads_no = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0
addprocs(threads_no)
@everywhere using Pkg
@everywhere Pkg.activate(".")


@everywhere using LinearAlgebra
@everywhere using BinaryOptimization
@everywhere using NPZ

@everywhere function simplified_tsp_qubo(w::Matrix{Float64})
    qubo = tsp_qubo(w, A=1.)
    n, _ = indexcapacities(vars(qubo))
    set_one!(qubo, [1,n])
    for i=1:(n-1)
        set_zero!(qubo, [1,i])
        set_zero!(qubo, [i+1,n])
    end
    qubo
end

@everywhere function simplified_tsp_hobo(w::Matrix{Float64})
    hobo = tsp_hobo(w, A=1.)
    n, logn = indexcapacities(vars(hobo))
    for (i, b) = enumerate(reverse(digits(n-1, base=2, pad=logn)))
        if b == 0
            set_zero!(hobo, [1,i])
        else
            set_one!(hobo, [1,i])
        end
    end

    if ispow2(n-1)
        for i=2:n
            set_zero!(hobo, [i, 1])
        end
    end
    hobo
end

@everywhere function hobo_emulated(qubo::Vector{<:Real}, n::Int)
    u = npzread("permutations_cutters/permcut_$(n-1).npy") # n-1 because of simplification
    true_u = reduce(kron, fill(u, n-1))
    true_u*qubo
end

println("arguments are fine. Generating data started ($threads_no threads, max $nmax cities)")


@sync @distributed for n = 3:nmax
    println(">>> $n ($nmax)")

    qubo = bom_to_hamiltonian(simplified_tsp_qubo(zeros(n, n)))
    hobo = bom_to_hamiltonian(simplified_tsp_hobo(zeros(n, n)))
    hobo_emu = hobo_emulated(qubo, n)

    npzwrite("$dir_out/qubo_hamilton_$n.npz", qubo)
    npzwrite("$dir_out/hobo_hamilton_$n.npz", hobo)
    npzwrite("$dir_out/hobo_emu_hamilton_$n.npz", hobo_emu)

    @assert length(hobo) == length(hobo_emu)
    @assert length(filter(iszero, hobo_emu)) == length(filter(iszero, hobo))
    @assert length(filter(iszero, hobo_emu)) == length(filter(iszero, qubo))

    if !(0. ≈ minimum(qubo) ≈ minimum(hobo) ≈ minimum(hobo_emu))
        println("Something wrong, minimal energy should be 0.")
    end
end

rmprocs()
