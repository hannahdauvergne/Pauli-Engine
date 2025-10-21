"""
List of function to compute several observables. 
 TODO: Comment it
"""
module Observables_Develop
export density, virial_energies, pair_correlation, 
       One_Body_Density_Matrix, One_Body_Density_Matrix_spatial,
       One_Body_Density_Matrix_eigvals, pair_correlation_spatial

include("Operators.jl")
include("Basis_creation.jl")
include("Hamiltonian_creation.jl")

using .Operators_Develop
using .BasisCreation_Develop
using .HamiltonianCreation_Develop
using SpecialPolynomials, Polynomials
using LinearAlgebra

function density(
    state :: Vector{Float64},
    particles :: Vector{Int},
    nho :: Int,
    type :: String,
    x :: Vector{Float64} = collect(-5.0:0.01:5.0)
    )

    basis, inv_base = do_basis(particles, nho, type)
    ndeg = size(particles, 1)
    nstates = size(basis, 1)
    OBDM = One_Body_Density_Matrix(state, particles, nho, type)
    density_function = zeros((ndeg+2,length(x)))
    density_function[1,:] = x
    for p in 1:ndeg
        for i in 1:nho
            for j in 1:nho
                density_function[p+1, :] += OBDM[p, i, j] * (HO_wavefunctions(i-1, x) .* HO_wavefunctions(j-1, x))
            end
        end
    end
    density_function[end, :] = sum(density_function[2:end-1, :], dims=1)[1, :]
    return density_function
end


function One_Body_Density_Matrix(
    state :: Vector{Float64},
    particles :: Vector{Int},
    nho :: Int,
    type :: String
    )

    basis, inv_basis = do_basis(particles, nho, type)
    ndeg = size(particles, 1)
    nstates = size(basis, 1)
    nsites = ndeg * nho
    
    if type== "fermi"
        op_an = fop_an
        op_cr = fop_cr
        maxE = nho + div(sum(particles .* (particles .+ 1)), 2) - maximum(particles)
    elseif type == "bose"
        op_an = bop_an
        op_cr = bop_cr
        maxE = nho
    else
        println("particletype must be fermi or bose")
        return
    end
    OBDM = zeros((ndeg,nho,nho))
    for i in 1:nstates
        state_0 = basis[i,:]
        for bi in unique(state_0)
            state_1, phase1 = op_an(copy(state_0), bi)
            E_2 = min(maxE - sum(div.(state_1 .- 1 , ndeg)),nho) #check
            for bk in (bi-1)%ndeg+1:ndeg:ndeg*E_2
                state_f,phase2 = op_cr(copy(state_1), bk)
                if haskey(inv_basis,state_f)
                    j = inv_basis[state_f]                        
                    phase = phase1*phase2
                    OBDM[bi%ndeg+1,div(bi-1,ndeg)+1,div(bk-1,ndeg)+1] += phase*state[i]*conj(state[j])
                end
            end
        end
    end
    return OBDM
end


function One_Body_Density_Matrix_spatial(
    state :: Vector{Float64},
    particles :: Vector{Int},
    nho :: Int,
    type :: String,
    x :: Vector{Float64} = collect(-5.0:0.01:5.0)
    )

    OBDM = One_Body_Density_Matrix(state, particles, nho, type)
    ndeg = size(particles, 1)
    npoints = length(x)
    OBDM_s = zeros((ndeg+2,npoints,npoints))
    for p in 1:ndeg
        for i in 1:nho
            for j in 1:nho
                OBDM_s[p+1, :, :] += OBDM[p, i, j] * (HO_wavefunctions(i-1, x) * HO_wavefunctions(j-1, x)')
            end
        end
    end
    OBDM_s[1,:,:] = x * x'
    OBDM_s[end, :, :] = sum(OBDM_s[2:end-1, :, :], dims=1)[1, :, :]
    return OBDM_s
end

function One_Body_Density_Matrix_eigvals(
    state :: Vector{Float64},
    particles :: Vector{Int},
    nho :: Int,
    type :: String
    )

    OBDM = One_Body_Density_Matrix(state, particles, nho, type)
    ndeg = size(particles, 1)
    vals = zeros(ndeg*nho)
    for i in 1:ndeg
        vals[(i-1)*nho+1:i*nho] = eigvals(OBDM[i,:,:])
    end
    sort!(vals,rev=true)
    return vals
end


function virial_energies(
    state :: Vector{Float64},
    particles :: Vector{Int},
    nho :: Int,
    type :: String,
    gint :: Vector{Float64}
    )
    #### Check if we can avoid create the Hamiltonians so many times
    H_kin = sum(one_body_Hamiltonian(particles, nho, type, etyp = "kin"))
    H_pot = sum(one_body_Hamiltonian(particles, nho, type, etyp="pot"))
    H_ID =  sum(gint .* two_body_Hamiltonian(particles, nho, type))

    E_kin = conj(state)' * H_kin * state
    E_pot = conj(state)' * H_pot * state
    E_int = conj(state)' * H_ID * state
    E_total = E_kin + E_pot + E_int
    virial = 2*E_kin - 2*E_pot + E_int

    return E_total,E_kin,E_pot,E_int,virial
end


function pair_correlation(
    state :: Vector{Float64},
    particles :: Vector{Int},
    nho :: Int,
    type :: String
    )

    basis, inv_basis = do_basis(particles, nho, type)
    ndeg = size(particles, 1)
    nstates = size(basis, 1)
    nsites = ndeg * nho

    if type== "fermi"
        op_an = fop_an
        op_cr = fop_cr
        maxE = nho + div(sum(particles .* (particles .+ 1)), 2) - maximum(particles)
    elseif type == "bose"
        op_an = bop_an
        op_cr = bop_cr
        maxE = nho
    else
        println("particletype must be fermi or bose")
        return
    end

    pair_corr = zeros((ndeg,ndeg,nho,nho, nho, nho))

    for i in 1:nstates
        state_0 = basis[i,:]
        for bi in unique(state_0)
            state_1, phase1 = op_an(copy(state_0), bi)
            for bj in unique(state_1[1:end-1])
                state_2, phase2 = op_an(copy(state_1), bj)
                E_2 = min(maxE - sum(div.(state_1 .- 1 , ndeg)),nho) #check
                for bk in (bj-1)%ndeg+1:ndeg:ndeg*E_2
                    state_3,phase3 = op_cr(copy(state_2), bk)
                    condition = min((E_2-div(bk,ndeg)+1) ,nho) #check (specially for fermions)
                    for bl in (bi-1)%ndeg+1:ndeg:condition*ndeg
                        state_f,phase4 = op_cr(copy(state_3), bl)
                        if haskey(inv_basis,state_f)
                            j = inv_basis[state_f]                        
                            phase = phase1*phase2*phase3*phase4
                            pair_corr[bi%ndeg+1,bj%ndeg+1,div(bi-1,ndeg)+1,div(bj-1,ndeg)+1,div(bk-1,ndeg)+1,div(bl-1,ndeg)+1] += phase*state[i]*conj(state[j])
                        end
                    end
                end
            end
        end
    end
    return pair_corr
end

# It is extremately slow....
function pair_correlation_spatial(
    state :: Vector{Float64},
    particles :: Vector{Int},
    nho :: Int,
    type :: String,
    x :: Vector{Float64} = collect(-5.0:0.01:5.0)
    )

    pair_corr = pair_correlation(state, particles, nho, type)
    ndeg = size(particles, 1)
    npoints = length(x)
    pair_corr_s = zeros((ndeg,ndeg,npoints,npoints))
    
    # Precompute all required HO wavefunctions
    ho_wfs = [HO_wavefunctions(n, x) for n in 0:nho-1]

    for p1 in 1:ndeg
        for p2 in 1:ndeg
            # Preallocate a temporary matrix for accumulation
            temp = zeros(npoints, npoints)
            for i in 1:nho
                ψi = ho_wfs[i]
                for j in 1:nho
                    ψj = ho_wfs[j]
                    for k in 1:nho
                        ψk = ho_wfs[k]
                        for l in 1:nho
                            ψl = ho_wfs[l]
                            coeff = pair_corr[p1, p2, i, j, k, l]
                            # Outer product of (ψi .* ψl) and (ψj .* ψk)
                            temp .+= coeff .* (ψi .* ψl) * ((ψj .* ψk)')
                        end
                    end
                end
            end
            pair_corr_s[p1, p2, :, :] = temp
        end
    end
        
    return pair_corr_s
end


function HO_wavefunctions(
    n :: Int,
    x :: Vector{Float64}
    )
    # ℏ = m = ω = 1 units
    norm = 1.0 / sqrt(sqrt(π) * 2.0^n * factorial(big(n)))
    ψ = norm * exp.(-x.^2 / 2) .* [basis(Hermite, n)(xi) for xi in x]
    return ψ
end


end
