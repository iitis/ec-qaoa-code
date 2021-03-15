using Distributed
using Pkg
Pkg.activate(".")

include("experiment_handler.jl")

# how many hobos/qubos
no_experiments = 40
# maximum k
max_k = 30
# how much each hamiltonian
repeating = 10

# how many procs
addprocs(threads_no)
@everywhere using Pkg
@everywhere Pkg.activate(".")
##
@everywhere using JLD2
@everywhere using FileIO
@everywhere using NPZ
@everywhere include("qaoa_optimizers.jl")

##

for _ in [1]
    upper = Dict("qubo" => 1. * pi, "hobo_emu" => 1. * pi, "hobo" => 2. * pi)
    function generate_experiment(mode::String, i::Int)
        println("######### experiment $i ($mode, $max_k) #########")
        kbreak = 5
        
        hamiltonian = npzread("$dir_in/$(mode)_$i.npz")
        
        results = Dict("experiment" => i,
                       "matrix_cost" => npzread("$dir_in/matrix_cost_$i.npz"),
                       "k" => max_k,
                       "mode" => mode,
                       "upper" => upper[mode])
        n = Int(log2(length(hamiltonian)))
        d = load_sparsers(n)
        for m = 1:repeating
            filename = "$dir_out/$mode-exp$i-m$m-result.jld2"
            if isfile(filename)
                println("file $filename exists")
                continue
            end
            results["small_results"] = qaoa(hamiltonian, kbreak-1, d=d, upper=upper[mode])
            t1 = time()
            results["large_results"] = qaoa_trajectories_periodic(hamiltonian, max_k, kmin=kbreak, d=d, upper=upper[mode])
            
            @show (mode, time() - t1)
            save(filename, results)
        end
        true
    end

    exp_tests = Iterators.product(["hobo", "qubo", "hobo_emu"], 1:no_experiments)

    pmap(i -> generate_experiment(i...), exp_tests)
end

rmprocs()
