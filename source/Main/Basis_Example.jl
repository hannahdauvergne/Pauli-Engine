include("Basis.jl")

function tiny_basis_creation(N::Int, lvls::Int, particletype::String)
    basis = Dict{Vector{Int}, Int}()
    index = 1

    # Loop over all possible ways of putting N particles into lvls levels
    for occs in Iterators.product(0:N, 0:N)
        if sum(occs) == N   # must have exactly N particles total
            occ = collect(occs)

            if particletype == "bose"
                basis[occ] = index
                index += 1
            elseif particletype == "fermi"
                # Fermions: only 0 or 1 per level
                if all(x -> x ≤ 1, occ)
                    basis[occ] = index
                    index += 1
                end
            else
                error("particletype must be 'bose' or 'fermi'")
            end
        end
    end

    return basis
end


basis_bose = tiny_basis_creation(2, 2, "bose")
for (state, idx) in basis_bose
    println("State $idx → $state")
end

