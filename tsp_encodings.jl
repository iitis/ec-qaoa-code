using BinaryOptimization

# hamiltonian cycle
function hamiltonian_cycle_qubo(n::Int, vars::Variables, A1::Real=1, A2::Real=A1)
    @assert vars.indexcapacities == [n, n]
    hobo = HOBO(vars)

    for t=1:n
        hobo_tmp = HOBO(vars)
        for j=1:n
            hobo_tmp += vars[t, j]
        end
        hobo += A1 * (1. - hobo_tmp)^2
    end

    for j=1:n
        hobo_tmp = HOBO(vars)
        for t=1:n
            hobo_tmp += vars[t, j]
        end
        hobo += A2 * (1. - hobo_tmp)^2
    end
    hobo
end

function hamiltonian_cycle_qubo(n::Int, name::String="Hamiltonian Cycle Comlete Graph")
    hamiltonian_cycle_qubo(n, Variables(name, [n, n], ["timestep", "vertex"]))
end

# tsp
function tsp_qubo(w::Matrix{T}, vars::Variables, B::Real=1., A1::Real=2*B*maximum(w), A2::Real=2*B*maximum(w)) where T<:Real
    n = size(w, 1)
    @assert size(w, 1) == size(w, 2)
    @assert vars.indexcapacities == [n, n]
    @assert A1 >= 0 && A2>=0 && B >= 0

    hobo_hamiltonian = hamiltonian_cycle_qubo(n, vars, A1, A2)
    hobo_tsp = HOBO(vars)
    #remove according to graph
    for u = 1:n, v = 1:n, t = 1:n
        if u != v
            hobo_tsp += w[u, v] * vars[t, u] * vars[t % n + 1, v]
        end
    end
    hobo_hamiltonian + B*hobo_tsp
end


function tsp_qubo(w::Matrix; name::String="TSP", B::Real=1., A::Real=2*B*maximum(w))
    tsp_qubo(w, Variables(name, collect(size(w)), ["timestep", "city"]), B, A, A)
end
