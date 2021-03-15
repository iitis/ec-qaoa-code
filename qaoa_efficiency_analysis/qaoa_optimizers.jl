using LinearAlgebra
using SparseArrays
using Optim
using Random
using NPZ
using LineSearches
include("sparse_generator_loader.jl")
include("qaoa_utils.jl")

import Optim: retract!, project_tangent!


function w_state(N::Int) # N - number of qubits
    result = zeros(ComplexF64, 2^N)
    # hack, works for small N only
    for k = [bitstring(2^i)[end-N+1:end] for i=0:(N-1)]
        ket_vec = [1]
        for b = k
            ket_vec = kron(ket_vec, b == '0' ? [1, 0] : [0, 1])
        end
        result += ket_vec
    end

    result/sqrt(N)
end

function w_init_state(t_dim, repeat)
    wstate = w_state(t_dim)
    result = [1]
    for _ = 1:repeat 
        result = kron(result, wstate)
    end
    result
end

struct Periodic <: Manifold
    periods::Vector{Float64}
    function Periodic(v::Vector{Float64})
        @assert all(v .> 0.)
        new(v)
    end
end

periods(p::Periodic) = p.periods

function retract!(p::Periodic, x)
    per = periods(p)
    neg_ind = findall(x -> x<zero(0.), x)
    correctors_int = ceil.(Int, abs.(x[neg_ind] ./ per[neg_ind]))
    x[neg_ind] .+= per[neg_ind] .* correctors_int
    x .%= per
    x
end

project_tangent!(p::Periodic, g, x) = g

##
function qaoa(hamiltonian::Vector{T},
              kmax::Int;
              upper::Float64=Float64(pi),
              d::Dict=load_sparsers(Int(log2(length(hamiltonian)))),
              init_state::Vector{ComplexF64}=fill(one(ComplexF64)/sqrt(length(hamiltonian)),length(hamiltonian))) where T<:Real
    @assert kmax >= 1
    @assert upper > 0.

    n = length(hamiltonian)
    tmp_data = Dict("state" => zeros(ComplexF64, n),
                    "mul_vec" => zeros(ComplexF64, n),
                    "v" => zeros(ComplexF64, 2*n),
                    "tmp_vec" => zeros(ComplexF64, n),
                    "tmp_vec2" => zeros(ComplexF64, n),
                    "d" => d)
    fg! = Optim.only_fg!((F, G, x) -> _fg!(F, G, hamiltonian, x, tmp_data, init_state))
    opt = Optim.Options(g_tol=1e-5, x_tol=1e-5, f_tol=1e-5, allow_f_increases=true, iterations=20_000)

    results = Any[]
    for k = 1:kmax
        push!(results, 0)
        periods = vcat(fill(upper, k), fill(pi, k))
        optimizer = LBFGS(manifold=Periodic(periods))#, alphaguess=InitialStatic(alpha=0.001))


        converged = false
        while !converged
            init_times = rand(2*k) .* periods
            results[end] = optimize(fg!, init_times, optimizer, opt)
            if Optim.converged(results[end])
                converged = true
            end
        end
    end
    results
end

function qaoa_trajectories_periodic(hamiltonian::Vector{T},
                                    kmax::Int;
                                    upper::Float64=2*pi,
                                    kmin::Int=5,
                                    d::Dict=load_sparsers(Int(log2(length(hamiltonian)))),
                                    level_repeat::Int=5,
                                    init_state::Vector{ComplexF64}=fill(one(ComplexF64)/sqrt(length(hamiltonian)),length(hamiltonian))) where T<:Real

    @assert kmax >= kmin >= 1
    @assert upper > 0.

    n = length(hamiltonian)
    tmp_data = Dict("state" => zeros(ComplexF64, n),
                    "mul_vec" => zeros(ComplexF64, n),
                    "v" => zeros(ComplexF64, 2*n),
                    "tmp_vec" => zeros(ComplexF64, n),
                    "tmp_vec2" => zeros(ComplexF64, n),
                    "d" => d)

    fg! = Optim.only_fg!((F, G, x) -> _fg!(F, G, hamiltonian, x, tmp_data, init_state))
    opt = Optim.Options(g_tol=1e-3, x_tol=1e-3, f_tol=1e-3, allow_f_increases=true, iterations=20_000)

    converged = false
    results = Any[]
    while !converged
        results = Any[]
        converged = true
        init_times = rand(2*kmin) .* vcat(fill(upper, kmin), fill(pi, kmin))
        best_t = copy(init_times) # only for initialization
        for k = kmin:kmax
            @show
            periods = vcat(fill(upper, k), fill(pi, k))
            optimizer = LBFGS(manifold=Periodic(periods))
            # optimizer = LBFGS(manifold=Periodic(periods), alphaguess=InitialStatic(alpha=0.001))            
            res_tmp = []
            
            for _ = 1:level_repeat
                if k > kmin
                    init_times = vcat(best_t[1:k], rand()*upper, best_t[k+1:end], rand()*pi)
                else
                    init_times = rand(2*kmin) .* vcat(fill(upper, kmin), fill(pi, kmin))
                end
                push!(res_tmp, optimize(fg!, init_times, optimizer, opt))
            end
            res_tmp = filter(Optim.converged, res_tmp)
            if length(res_tmp) == 0
                println("Failed qaoa_trajectories_periodic $k $level_repeat")
                converged = false
                break
            end
            _, pos = findmin(Optim.minimum.(res_tmp))
            push!(results, res_tmp[pos])
            best_t = Optim.minimizer(results[end])
        end

    end

    results
end
