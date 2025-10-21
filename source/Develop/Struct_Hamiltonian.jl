include("Operators.jl")
include("Basis_creation.jl")
include("MintCreation.jl")


using .Operators_Develop
using .BasisCreation_Develop
using .MintCreation_Develop

""" 
    Calculation of a one body Hamiltonian matrix. It supports three options:
        "ho" for the single-particle harmonic oscillator Hamiltonian.
        "kin" for the single-particle kinetic part.
        "pot" fot the single-particle harmonic oscillator part.
    It modifies the internal stored Hamiltonians o the object.
    Is possible to save and load computed Hamiltonians.
    It creates an array of matrices, one per internal component.
"""
function one_body_Hamiltonian!(
    system :: Few_Particles_Hamonic_Oscillator; 
    etyp :: String = "ho"
    )
    particles = system.system.particles
    sites = system.system.lvls
    ptype = system.system.particletype

    if system.save_local 
        if etyp == "ho"
            file = system.path * system.path_H * "Hho_" * system.system.id * ".jl2d"
            if isfile(file)
                @load file H
                system.H_mat =Hamiltonians(H,system.H_mat.HI,system.H_mat.Hkin,system.H_mat.Hpot)
                return
            end
        elseif etyp == "kin"
            file = system.path * system.path_H * "Hkin_" * system.system.id * ".jl2d"
            if isfile(file)
                @load file H
                system.H_mat =Hamiltonians(system.H_mat.H0,system.H_mat.HI,H,system.H_mat.Hpot)
                return
            end
        elseif etyp == "pot"
            file = system.path * system.path_H * "Hpot_" * system.system.id * ".jl2d"
            if isfile(file)
                @load file H
                system.H_mat =Hamiltonians(system.H_mat.H0,system.H_mat.HI,system.H_mat.Hkin,H)
                return
            end
        end
    end


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
    basis = system.system.basis
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
                    push!(value[bi%deg+1],phase*0.25*sqrt(ival*(ival-1)))
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
                    push!(value[bi%deg+1],phase*0.25*sqrt((ival+2)*(ival+1)))
                end
            end
        end
    end
    H = [sparse(i_index[i],j_index[i],value[i],dim, dim) for i in 1:deg]
    if system.save_local
        @save file H
    end

    if etyp == "ho"
        system.H_mat =Hamiltonians(H,system.H_mat.HI,system.H_mat.Hkin,system.H_mat.Hpot)
    elseif etyp == "kin"
        system.H_mat =Hamiltonians(system.H_mat.H0,system.H_mat.HI,H,system.H_mat.Hpot)
    elseif etyp == "pot"
        system.H_mat =Hamiltonians(system.H_mat.H0,system.H_mat.HI,system.H_mat.Hkin,H)
    end
end

"""
    Creates a two particle Hamiltonian matrix. For now there is only one possibility.
    "cont" creates a delta contact interaction between the particles. 
    It creates an array of matrices, each one corresponding to a pair of the different pair of components.
    It is sorted as: [H_{1,1}, H_{1,2}, ... , H_{1,N}, H_{2,2}, H_{2,3}, ... , H_{2,N}, H_{3,3}, ... ,  H_{3,N}, ... , H_{N,N} ],
        where H_{i,j} is the interacting Hamiltonian between the components i and j (note that H_{i,j} = H_{j,i}) and N is the numbe of components.
"""
function two_body_Hamiltonian!(
    system :: Few_Particles_Hamonic_Oscillator ;
    vtype :: String = "cont"
    )
    particles = system.system.particles
    sites = system.system.lvls
    ptype = system.system.particletype

    if system.save_local
    file = system.path * system.path_H * "Hint_" * system.system.id * ".jl2d"
        if isfile(file)
            @load file H
            system.H_mat =Hamiltonians(system.H_mat.H0,H,system.H_mat.Hkin,system.H_mat.Hpot)
            return
        end
    end

    if ptype == "fermi"
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

    basis = system.system.basis
    spin = length(particles)
    deg = Int(spin*(spin+1)/2) #Is different for bosons and fermions. This is for bososn, and then some components will be zero
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
    if system.save_local
        @save file H
    end

    system.H_mat =Hamiltonians(system.H_mat.H0,H,system.H_mat.Hkin,system.H_mat.Hpot)
end

"""
    Transforms the single-particle and the two-particle Hamiltonians to the total Hamiltonian for a specific itneraction strength
"""
function total_Hamiltonian(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    g :: Union{Vector{Float64},Float64}
    )

    ideg_int = size(Hi)[1]

    if length(g) != ideg_int
        g = g[1] * ones(ideg_int)
    end

    H= sum(H0)+sum(g .*Hi)
    return H
end


