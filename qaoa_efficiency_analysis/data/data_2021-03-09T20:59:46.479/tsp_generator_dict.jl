using Pkg
Pkg.activate(".")



if length(ARGS) == 0 || ARGS[1] == "--help" || ARGS[1] == "-h" || ARGS[1] == "help"
    println("\n\tFormat: julia tsp_generator_dict.jl  n number

n is the number of cities, number is number of instances.
0th instance is Hamilton problem.
n should be small. Data are saved to 'tsp_dict'.
")
    exit()
end

n = parse(Int, ARGS[1])
@assert n >= 1

include("qaoa_utils.jl")
include("sparse_generator_loader.jl")
if !isdir("sparse_1qubitgate_data")
    mkdir("sparse_1qubitgate_data")
end
#generator(n^2)

repeating = parse(Int, ARGS[2])
@assert repeating >= 1

using BinaryOptimization
using LinearAlgebra
using NPZ

function save_dict(d::Dict, filename::String)
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

function simplified_tsp_qubo(w::Matrix{Float64}, A::Float64=maximum(w))
    qubo = tsp_qubo(w,A=A)
    n, _ = indexcapacities(vars(qubo))
    set_one!(qubo, [n,n])
    for i=1:(n-1)
        set_zero!(qubo, [n,i])
        set_zero!(qubo, [i,n])
    end
    qubo
end

println("arguments are fine. Generating data started ($repeating samples)")
for m = 1:repeating
    println("> $m ($repeating)")
    w = rand(1:10, (n,n))
    w += Matrix{Float64}(transpose(w) - 2*diagm(diag(w)))

    qubo = Ising(simplified_tsp_qubo(w))
    #qubo_vec = bom_to_hamiltonian(qubo)
    #vec_n = length(qubo_vec)
    #qubits_no = n^2

    save_dict(qubo.values, "tsp_dict/tsp_$n-$m.npz")
    npzwrite("tsp_dict/matrix_cost_$n-$m.npz", w)
    #npzwrite("tsp_dict/qubo_vec_$n-$m.npz", qubo_vec)

    #=
    angles = rand(10)*2*pi

    tmp_data = Dict("state" => zeros(ComplexF64, vec_n),
                    "mul_vec" => zeros(ComplexF64, vec_n),
                    "v" => zeros(ComplexF64, 2*vec_n),
                    "tmp_vec" => zeros(ComplexF64, vec_n),
                    "tmp_vec2" => zeros(ComplexF64, vec_n),
                    "d" => load_sparsers(qubits_no))

    state = fill(ComplexF64(1.) / sqrt(vec_n), vec_n)
    state = _state!(qubo_vec, angles, tmp_data, state)
    npzwrite("tsp_dict/angles_test_$n-$m.npz", angles)
    npzwrite("tsp_dict/state_test_$n-$m.npz", state)
    =#
end
