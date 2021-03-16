using Pkg
Pkg.activate(".")
using Distributed
try
    addprocs(parse(Int, ARGS[1]))
catch
    nothing
end
##
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using PyPlot
@everywhere using JLD2, FileIO
@everywhere using LinearAlgebra
@everywhere using SparseArrays
@everywhere using Optim
@everywhere using Statistics
@everywhere using NPZ
@everywhere include("sparse_generator_loader.jl")
@everywhere include("qaoa_optimizers.jl")


rc("text", usetex=true)
rc("font", family="serif", size=8)

cities_no = [4]
data_dir = Dict(4 => "data/data_2021-03-09T20:59:46.479/data/")
data_instances = Dict(4 => "data/data_2021-03-09T20:59:46.479/tsp4/")
k_max = 30
k_min = 1
experiment = 40 
repeating = 10
## lowe eigenstate probability

function reverse_hamiltonian_thr(v::Vector{Float64}, thr::Float64)
    map(x -> x == 0 ? 1. : 0., v)
end

function up_low_std(data::Vector{Float64})
    the_mean = mean(data)
    up_std = filter(x-> x > the_mean, data)
    up_std = length(up_std) <= 2 ? 0 : sqrt(sum((up_std .- the_mean).^2)/length(up_std))
    low_std = filter(x-> x < the_mean, data)
    low_std = length(low_std) <= 2 ? 0 : sqrt(sum((low_std .- the_mean).^2)/length(low_std))
    the_mean, up_std, low_std
end

# hobo_6_result_14.jld2
function get_best(filename::String)
    d = load(filename)["results"]
    output = d[findmin(minimum.(d))[2]]
end

##

if "generate" ∈ ARGS
    for city_no = cities_no
        println("######## city $city_no ##########")
        for mode = ["qubo", "hobo", "hobo_emu"]
            ham_hamilton = npzread("hamilton/$(mode == "hobo_emu" ? mode : mode[1:4])_hamilton_$city_no.npz")
            ham_feasible = reverse_hamiltonian_thr(ham_hamilton, 0.01)
            n = length(ham_hamilton)
            qubits_no = Int(log2(n))

            function generate_data(exp::Int)
                println("$city_no $mode $exp")
                tmp_data = Dict("state" => zeros(ComplexF64, n),
                                "mul_vec" => zeros(ComplexF64, n),
                                "v" => zeros(ComplexF64, 2*n),
                                "tmp_vec" => zeros(ComplexF64, n),
                                "tmp_vec2" => zeros(ComplexF64, n),
                                "d" => load_sparsers(qubits_no))
                ham_tsp = npzread("$(data_instances[city_no])/$(mode)_$exp.npz")
                ham_tsp_feasible = ham_tsp .* ham_feasible
                true_min = minimum(collect(filter(x->x>0, ham_tsp_feasible)))
                true_max = maximum(collect(filter(x->x>0, ham_tsp_feasible)))

                data_prob = zeros(k_max-k_min+1, repeating)
                data_energy = zeros(k_max-k_min+1, repeating)
                for m = 1:repeating
                    sol = load("$(data_dir[city_no])/$(mode)-exp$exp-m$m-result.jld2")
                    results = vcat(sol["small_results"], sol["large_results"])
                    @assert length(results) >= k_max-k_min+1 
                    @assert all(Optim.converged.(results))
                    for (ind, el) = enumerate(results)
                       
                        state = fill(ComplexF64(1.) / sqrt(n), n)
                        data_prob[ind, m] = _energy_diffham!(ham_feasible, ham_tsp, Optim.minimizer(el), tmp_data, state)

                        state = fill(ComplexF64(1.) / sqrt(n), n)
                        state = _state!(ham_tsp, Optim.minimizer(el), tmp_data, state)
                        state .= state .* ham_feasible
                        state .=  state ./ norm(state)
                        energy_state = sum(abs2.(state) .* ham_tsp_feasible)
                        data_energy[ind, m] = (energy_state - true_min)/(true_max - true_min)

                    end
                end
                npzwrite("$(data_dir[city_no])/$mode-city$city_no-exp$exp-plotdata-prob.npz", data_prob)
                npzwrite("$(data_dir[city_no])/$mode-city$city_no-exp$exp-plotdata-energy.npz", data_energy)
                return true
            end
            pmap(generate_data, 1:experiment)
        end
    end
end

if "plot" ∈ ARGS
    # probs
    city_no = 4

    fig, ax = subplots(figsize=[2.2,1.4], nrows=1, ncols=1, sharex=true, sharey=true)
    println("######## city $city_no ##########")
    data_x = k_min:k_max
    for (color, mode, label) = zip(["r:", "b--", "k-"], ["hobo", "hobo_emu", "qubo"], ["HOBO", "enc-change", "QUBO-QAOA"])
        data_y = zeros(length(data_x), 0)
        for exp = 1:experiment
            data_prob = npzread("$(data_dir[city_no])/$mode-city$city_no-exp$exp-plotdata-prob.npz")
            data_y = hcat(data_y, mapslices(maximum, data_prob, dims=[2]))
        end
        mean_data = mapslices(mean, data_y, dims=[2])[:]
        up_border = mean_data .+ mapslices(std, data_y, dims=[2])[:] 
        bot_border = mean_data .- mapslices(std, data_y, dims=[2])[:]  

        ax.plot(data_x, mean_data, "$color", label=label,markersize=.5)
        ax.fill_between(data_x, up_border, bot_border, color=color[[1]], alpha = 0.2)

        ax.vlines(5, 0., 1., color="k", linestyle="--", linewidth=.5)
        ax.set_xlabel("Number of levels")
    end
    ax.set_ylim(-0.05, 1.05)
    #ax.legend(loc=1, prop=Dict("size" => 6), bbox_to_anchor=[1.52,1.04])
    ax.set_ylabel("prob. of measuring\n in feasible space")  
    savefig("plots/feasible_prob_tsp_trajectories_multi.pdf", bbox_inches="tight")
    cla()

    # energy
    fig, ax = subplots(figsize=[2.2,1.4], nrows=1, ncols=1, sharex=true, sharey=true)
        
    println("######## city $city_no ##########")
    data_x = k_min:k_max
    for (color, mode, label) = zip(["r:", "b--", "k-"], ["hobo", "hobo_emu", "qubo"], ["HOBO-QAOA", "EC-QAOA", "QUBO-QAOA"])
        data_y = zeros(length(data_x), 0)
        for exp = 1:experiment
            data_prob = npzread("$(data_dir[city_no])/$mode-city$city_no-exp$exp-plotdata-energy.npz")
            data_y = hcat(data_y, mapslices(maximum, data_prob, dims=[2]))
        end
        mean_data = mapslices(mean, data_y, dims=[2])[:]
        up_border = mean_data .+ mapslices(std, data_y, dims=[2])[:] 
        bot_border = mean_data .- mapslices(std, data_y, dims=[2])[:]  

        ax.plot(data_x, mean_data, "$color", label=label, markersize=1, linewidth=1)
        ax.fill_between(data_x, up_border, bot_border, color=color[[1]], alpha = 0.2)
        ax.vlines(5, 0., 1., color="k", linestyle="--", linewidth=.5)
    end
    ax.legend(loc=1, prop=Dict("size" => 6), bbox_to_anchor=[1.62,1.03])
    ax.set_ylim(-0.05, 1.05)
    setp(ax, ylabel="rescaled energy\nfor feasible space")

    ax.set_xlabel("Number of levels")
    savefig("plots/feasible_energy_tsp_trajectories_multi.pdf", bbox_inches="tight")
    cla()
end