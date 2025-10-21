""" 
    Computes the density profile for a given state and interaction strenght.
        it computes the eigenstate requiered (if not previous computed) or load it.
        By computing (or loaing) the One BodyDensity Matrix, it computes the density profile for a grid x 
"""
function density_profile(
    system :: Few_Particles_Hamonic_Oscillator;
    x :: Vector{Float64} = collect(-5.0:0.01:5.0)
    )

    params = (system.target_state, system.target_g)
    if params != system.params.prev_params
        eigenstate!(system)
    end

    particles = system.system.particles
    nho = system.system.lvls

    ndeg = size(particles, 1)
    OBDM = One_Body_Density_Matrix(system) 
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
    dens = density_struct(x,density_function[2:end-1,:],density_function[end,:],system.target_g)
    return dens
end

"""
    Computes the One Body Density Matrix fora given eigenstate. 
"""
function One_Body_Density_Matrix(
    system :: Few_Particles_Hamonic_Oscillator
    )

    params = (system.target_state, system.target_g)
    sp = system.params
    if params != system.params.prev_params
        system.params = Internal_params(sp.prev_params,sp.target_eigenstate,sp.target_gcalc,Array{Float64,3}(undef, 0, 0, 0))
        eigenstate!(system)
    end

    if length(system.params.OBDM) != 0
        return system.params.OBDM
    end

    state = system.params.target_eigenstate 
    particles = system.system.particles
    nho = system.system.lvls
    type = system.system.particletype

    basis = system.system.basis
    ndeg = size(particles, 1)
    
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
    for state_0 in keys(basis)
        i = basis[state_0]
        for bi in unique(state_0)
            state_1, phase1 = op_an(copy(state_0), bi)
            E_2 = min(maxE - sum(div.(state_1 .- 1 , ndeg)),nho) #check
            for bk in (bi-1)%ndeg+1:ndeg:ndeg*E_2
                state_f,phase2 = op_cr(copy(state_1), bk)
                if haskey(basis,state_f)
                    j = basis[state_f]                        
                    phase = phase1*phase2
                    OBDM[(bi-1)%ndeg+1,div(bi-1,ndeg)+1,div(bk-1,ndeg)+1] += phase*state[i]*conj(state[j])
                end
            end
        end
    end
    system.params = Internal_params(sp.prev_params,sp.target_eigenstate,sp.target_gcalc,OBDM)
    return OBDM
end

""" 
    Compute the OBDM(x,x') for each component of the system
"""
function One_Body_Density_Matrix_spatial(
    system :: Few_Particles_Hamonic_Oscillator;
    x :: Vector{Float64} = collect(-5.0:0.01:5.0)
    )

    params = (system.target_state, system.target_g)
    if params != system.params.prev_params
        eigenstate!(system)
    end

    particles = system.system.particles
    nho = system.system.lvls

    OBDM = One_Body_Density_Matrix(system)
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

"""
    Compute the OBDM eigenvalues for each component, and then mix all of them and sort it
"""
function One_Body_Density_Matrix_eigvals(
    system :: Few_Particles_Hamonic_Oscillator
    )
    
    params = (system.target_state, system.target_g)
    if params != system.params.prev_params 
        eigenstate!(system)
    end

    particles = system.system.particles
    nho = system.system.lvls

    OBDM = One_Body_Density_Matrix(system)
    ndeg = size(particles, 1)
    vals = zeros(ndeg*nho)
    for i in 1:ndeg
        vals[(i-1)*nho+1:i*nho] = eigvals(OBDM[i,:,:])
    end
    sort!(vals,rev=true)
    return vals
end

""" 
    Computes the different contributions of the energy: Kinetic, potential and interacting energies.
    Gives the discrepancy with the virial theorem ( 2T + E_i = 2V ) as 2T - 2V + E_i.
    And the total energy computed as the sum of the components, to check with the diagonalization value (it should be the same).
"""
function virial_energies(
    system :: Few_Particles_Hamonic_Oscillator
    )
    #### Check if we can avoid create the Hamiltonians so many times

    if length(system.H_mat.Hkin) == 0
        one_body_Hamiltonian!(system, etyp = "kin")
    end
    if length(system.H_mat.Hpot) == 0
        one_body_Hamiltonian!(system,etyp = "pot")
    end
    if length(system.H_mat.HI) == 0
        two_body_Hamiltonian!(system)
    end
    params = (system.target_state, system.target_g)
    if params != system.params.prev_params
        eigenstate!(system)
    end
    
    state = system.params.target_eigenstate 
    H_kin = sum(system.H_mat.Hkin)
    H_pot = sum(system.H_mat.Hpot)
    H_ID =  sum(system.params.target_gcalc .* system.H_mat.HI) 

    E_kin = conj(state)' * H_kin * state
    E_pot = conj(state)' * H_pot * state
    E_int = conj(state)' * H_ID * state
    E_total = E_kin + E_pot + E_int
    virial = 2*E_kin - 2*E_pot + E_int

    return E_total,E_kin,E_pot,E_int,virial
end

"""
    Wavefunction of the Harmonic oscillator functions 
"""
function HO_wavefunctions(
    n :: Int,
    x :: Vector{Float64}
    )
    # ℏ = m = ω = 1 units
    norm = 1.0 / sqrt(sqrt(π) * 2.0^n * factorial(big(n)))
    ψ = norm * exp.(-x.^2 / 2) .* [basis(Hermite, n)(xi) for xi in x]
    return ψ
end