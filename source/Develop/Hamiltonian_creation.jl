
################
### Implement only the dictionary method. Seems to need exactly the same time
###############

################ Compared with the python code (In fact, this julia code is more optimized) 
################              for the case of particles = [1,1,1] and sites = 30
################                         in python needs:   92.26 seconds
################                          in julia needs:    3.44 seconds
################              for the case of particles = [1,1,1] and sites = 45
################                         in python needs:    600 seconds 
################                          in julia needs:     24 seconds
################  For the moment, the Backup code is even slower than the python version (and compute it as dense matrices)
module HamiltonianCreation_Develop
export one_body_Hamiltonian, two_body_Hamiltonian
include("Operators.jl")
include("Basis_creation.jl")
include("MintCreation.jl")
using SpecialFunctions
using SparseArrays

using .Operators_Develop
using .BasisCreation_Develop
using .MintCreation_Develop

"""
    one_body_Hamiltonian(particles, sites, ptype, etyp) => H
Computes a one body Hamiltonian for a certain parameters, as particles and harmonic oscillator modes.
It is implemented for different Hamiltonians.

# Inputs

    particles :: Vector{Int}        Number of particles in the system
    sites :: Int                    Number of modes used in the basis
    ptype :: String                 Type of particles, "bose" for bosons and "fermi" for fermions.

# Optional Inputs

    etyp :: String = "ho"           Hamiltonian created. Default the single-particle Harmonic Oscillator Hamiltonian ("ho"). 
                                    Also are implemented the kinetic ("kin") and the potential ("pot") independently.

# Returns

    H :: Vector{SparseMatrixCSC{Float64, Int}}          One-body Hamiltonian. The first index correspond to a different internal 
                                                        state of the particles.

"""
function one_body_Hamiltonian(
    particles :: Vector{Int} , 
    sites :: Int, 
    ptype :: String; 
    etyp :: String = "ho"
    )

    if ptype == "fermi"
        op_an = fop_an
        op_cr = fop_cr
    elseif ptype == "bose"
        op_an = bop_an
        op_cr = bop_cr
    else
        printl("particletype must be fermi or bose")
        return
    end
    basis = do_basis(particles, sites, ptype)
    deg = length(particles)
    dim = length(keys(basis))

    i_index = [ (Int[]) for _ in 1:deg ]
    j_index = [ (Int[]) for _ in 1:deg ]
    value = [ (Float64[]) for _ in 1:deg ]
    if etyp == "ho"
        for state in keys(basis)
            i = basis[state]
            for k in 1:deg
                push!(i_index[k],i)
                push!(j_index[k],i)
                push!(value[k],sum(div(val-1, deg) + 0.5 for val in state if (val-1) % deg == k-1))
            end
        end
    elseif etyp == "kin"
        for state_0 in keys(basis)
            i = basis[state_0] 
            for bi in unique(state_0)
                ival = div(bi-1,deg)
                state_1, phase1 = op_an(copy(state_0), bi)
                state_f,phase2 = op_cr(copy(state_1), bi-2*deg)
                if haskey(basis,state_f)
                    j = basis[state_f]                        
                    phase = phase1*phase2
                    push!(i_index[bi%deg+1],i)
                    push!(j_index[bi%deg+1],j)
                    push!(value[bi%deg+1],-phase*0.25*sqrt(ival*(ival-1)))
                end
                state_f,phase2 = op_cr(copy(state_1), bi)
                if haskey(basis,state_f)
                    j = basis[state_f]                        
                    phase = phase1*phase2
                    push!(i_index[bi%deg+1],i)
                    push!(j_index[bi%deg+1],j)
                    push!(value[bi%deg+1],phase*0.25*(2*ival +1))
                end
                state_f,phase2 = op_cr(copy(state_1), bi+2*deg)
                if haskey(basis,state_f)
                    j = basis[state_f]                        
                    phase = phase1*phase2
                    push!(i_index[bi%deg+1],i)
                    push!(j_index[bi%deg+1],j)
                    push!(value[bi%deg+1],-phase*0.25*sqrt((ival+2)*(ival+1)))
                end
            end
        end
    elseif etyp == "pot"
        for state_0 in keys(basis)
            i = basis[state_0]
            for bi in unique(state_0)
                ival = div(bi-1,deg)
                state_1, phase1 = op_an(copy(state_0), bi)
                state_f,phase2 = op_cr(copy(state_1), bi-2*deg)
                if haskey(basis,state_f)
                    j = basis[state_f]                        
                    phase = phase1*phase2
                    push!(i_index[bi%deg+1],i)
                    push!(j_index[bi%deg+1],j)
                    push!(values[bi%deg+1],phase*0.25*sqrt(ival*(ival-1)))
                end
                state_f,phase2 = op_cr(copy(state_1), bi)
                if haskey(basis,state_f)
                    j = basis[state_f]                        
                    phase = phase1*phase2
                    push!(i_index[bi%deg+1],i)
                    push!(j_index[bi%deg+1],j)
                    push!(values[bi%deg+1],phase*0.25*(2*ival +1))
                end
                state_f,phase2 = op_cr(copy(state_1), bi+2*deg)
                if haskey(basis,state_f)
                    j = basis[state_f]                        
                    phase = phase1*phase2
                    push!(i_index[bi%deg+1],i)
                    push!(j_index[bi%deg+1],j)
                    push!(values[bi%deg+1],phase*0.25*sqrt((ival+2)*(ival+1)))
                end
            end
        end
    end
    H = [sparse(i_index[i],j_index[i],value[i],dim, dim) for i in 1:deg]
    return H
end

"""
    two_body_Hamiltonian(particles, sites, ptype, vtype) => H
Computes a two body Hamiltonian for a certain parameters, as particles and harmonic oscillator modes.
It is implemented for different Hamiltonians.

# Inputs

    particles :: Vector{Int}        Number of particles in the system
    sites :: Int                    Number of modes used in the basis
    ptype :: String                 Type of particles, "bose" for bosons and "fermi" for fermions.

# Optional Inputs

    vtype :: String = "cont"         Type of Hamiltonian created. Default the contact interaction ("cont")

# Returns

    H :: Vector{SparseMatrixCSC{Float64, Int}}          Two-body Hamiltonian. The first index correspond to a different 
                                                        pair of internal states of the particles. (need clarification)
"""
function two_body_Hamiltonian(
    particles :: Vector{Int} ,
    sites :: Int ,
    ptype :: String ;
    vtype :: String = "cont"
    )
    if ptype== "fermi"
        op_an = fop_an
        op_cr = fop_cr
        maxE = sites + div(sum(particles .* (particles .+ 1)), 2) - maximum(particles) 
    elseif ptype == "bose"
        op_an = bop_an
        op_cr = bop_cr
        maxE = sites 
    else
        println("particletype must be fermi or bose")
        return
    end

    mint = compute_mint_recursive_self_v2(sites)

    basis = do_basis(particles,sites,ptype)
    spin = length(particles)
    deg = Int(spin*(spin+1)/2)
    dim = length(keys(basis))
    i_index = [ (Int[]) for _ in 1:deg ]
    j_index = [ (Int[]) for _ in 1:deg ]
    value = [ (Float64[]) for _ in 1:deg ]
    if vtype == "cont"
        for state in keys(basis)
            i = basis[state]
            for bi in unique(state)
                state_1, phase1 = op_an(copy(state), bi)
                for bj in unique(state_1[1:end-1])
                    state_2,phase2 = op_an(copy(state_1), bj)
                    E_2 = min(maxE - sum(div.(state_2 .- 1 , spin)),sites) 
                    for bk in (bj-1)%spin+1:spin:spin*E_2 
                        state_3,phase3 = op_cr(copy(state_2), bk)
                        condition = min((E_2-div(bk,spin)+1) ,sites) 
                        pair = ((div((bi-1),spin) + div((bj-1),spin) + div((bk-1),spin) ) % 2) * spin
                        for bl in (bi-1)%spin+1+pair:2*spin:condition*spin 
                            state_f,phase4 = op_cr(copy(state_3), bl)
                            if haskey(basis,state_f)
                                j = basis[state_f]                        
                                phase = phase1*phase2*phase3*phase4
                                
                                el = div.([bi,bj,bk,bl] .- 1,spin) .+1
                                sort!(el,rev=true)
                                lab = sort([(bi-1)%spin,(bj-1)%spin])
                                degcond = Int(lab[2]-lab[1]*(1+lab[1]-2*spin)/2)+1

                                push!(i_index[degcond],i)
                                push!(j_index[degcond],j)
                                push!(value[degcond],1/2*phase*mint[el[1],el[2],el[3],el[4]])
                            end
                        end
                    end
                end
            end
        end
    end
    H = [sparse(i_index[i],j_index[i],value[i],dim, dim) for i in 1:deg]
    return H 
end

function total_Hamiltonian(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    g :: Union{Vector{Float64},Float64}
    )

    ideg_int = size(Hi)

    if length(g) != ideg_int[1]
        g = g[1] * ones(ideg_int[1])
    end

    H= sum(H0)+sum(g .*Hi)
    return H
end

end
