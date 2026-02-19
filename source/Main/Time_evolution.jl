
function quench(
    system :: Few_Particles_Hamonic_Oscillator,
    g0 :: Float64,
    gf :: Float64
    )
    if length(system.H_mat.H0) == 0
        one_body_Hamiltonian!(system)
    end

    if length(system.H_mat.HI) == 0
        two_body_Hamiltonian!(system)
    end
    H0 = system.H_mat.H0
    Hi = system.H_mat.HI

    sizeH = size(H0[1])[1]
    v0 = rand(sizeH)

    H = total_Hamiltonian(H0,Hi,g0)
    decomp, history = partialschur(
        H,
        v1 = v0,
        nev = 5,
        which = :SR,
        tol = 1e-15,
    )
    # The system uses the supoerposition of the eigenstates of the prevous step to initialize the following calculation to reduce the computation time.
    v0 = vec(sum(decomp.Q; dims=2))
    state_0 = Matrix(decomp.Q)[:,1]

    overlap = 0.0
    numstates = 20
    projections = 0
    while (1 - overlap) > 5e-3
        H = total_Hamiltonian(H0,Hi,gf)
        if numstates > length(system.system.basis)
            println("maximum reached, need a larger basis")
            numstates = length(system.system.basis)
        end
        decomp, history = partialschur(
            H,
            v1 = v0,
            nev = numstates,
            which = :SR,
            tol = 1e-15,
        )
        projections = Matrix(decomp.Q)' * state_0
        overlap = sum(abs2.(projections))
        numstates += 20
    end
    indxs = abs2.(projections) .> 1e-5
    return (projections[indxs],real.(decomp.eigenvalues)[indxs],Matrix(decomp.Q)[:,indxs])
end

function quench_continuos(
    system :: Few_Particles_Hamonic_Oscillator,
    state_0,
    gf :: Float64;
    type = "int",
    omegaf :: Float64
    )
    
    H0 = system.H_mat.H0
    Hi = system.H_mat.HI
    Hb = system.H_mat.Hpot

    sizeH = size(H0[1])[1]
    v0 = rand(sizeH)

    overlap = 0.0
    numstates = 20
    projections = 0
    decomp = 0
    while (1 - overlap) > 5e-3
        if type == "int"
            H = total_Hamiltonian(H0,Hi,gf)
        elseif type == "breath"
            Ht = total_Hamiltonian(H0,Hi,gf)
            H = total_Hamiltonian([Ht],Hb,omegaf)
        end
        if numstates > length(system.system.basis)
            println("maximum reached, need a larger basis")
            numstates = length(system.system.basis)
        end
        decomp, history = partialschur(
            H,
            v1 = v0,
            nev = numstates,
            which = :SR,
            tol = 1e-15,
        )
        projections = Matrix(decomp.Q)' * state_0
        overlap = sum(abs2.(projections))
        numstates += 20
    end
    indxs = abs2.(projections) .> 1e-5
    return (projections[indxs],real.(decomp.eigenvalues)[indxs],Matrix(decomp.Q)[:,indxs])
end

function time_evolution(
    system :: Few_Particles_Hamonic_Oscillator,
    gvalues,
    times;
    type = "int",
    omegas = zeros(1)
    )
    if length(system.H_mat.H0) == 0
        one_body_Hamiltonian!(system)
    end

    if length(system.H_mat.HI) == 0
        two_body_Hamiltonian!(system)
    end

    steps = size(times)[1]
    xs = collect(range(-5,5,1001))
    Energy = zeros(steps)
    dens = zeros((steps,1001))
    norm = zeros(steps)
    fidelity = zeros(steps)

    #i=2
    #dt = times[i]-times[i-1]
    #proj,energies,states = quench(system,gvalues[i-1],gvalues[i])
    #exp_matrix = exp.(-1im .* energies .* dt)
    #tvec = exp_matrix' .* proj'
    #state = vec(tvec * states')

    #Energy[1] = energies[1]

    H0 = system.H_mat.H0
    if type == "int"
        Hi = system.H_mat.HI
    elseif type == "breath"
        Hi = system.H_mat.HI
        if length(system.H_mat.Hpot) == 0
            one_body_Hamiltonian!(system,etyp="pot")
        end
        Hb = system.H_mat.Hpot
    end
    if type == "int"
        H = total_Hamiltonian(H0,Hi,gvalues[1])
    elseif type == "breath"
        Ht = total_Hamiltonian(H0,Hi,gvalues[1])
        H = total_Hamiltonian([Ht],Hb,omegas[1])
    end

    sizeH = size(H0[1])[1]
    v0 = rand(sizeH)

    decomp, history = partialschur(
        H,
        v1 = v0,
        nev = 5,
        which = :SR,
        tol = 1e-15,
    )
    state = Matrix(decomp.Q)[:,1]


    Energy[1] = real.(state' * H * state)
    dens[1,:] = density_profile_time(system,state,x=xs)
    norm[1] = sum(abs2.(state))
    state_init = copy(state)
    fidelity[1] = abs.(state' * state_init)

    for i in 2:steps
        dt = times[i]-times[i-1]
        if type == "int"
            proj,energies,states = quench_continuos(system,state,gvalues[i])
        elseif type == "breath"
            proj,energies,states = quench_continuos(system,state,gvalues[i],type="breath",omegaf=omegas[i])
        end
        exp_matrix = exp.(-1im .* energies .* dt)
        tvec = exp_matrix' .* proj'
        state = vec(tvec * states')
        if type == "int"
            H = total_Hamiltonian(H0,Hi,gvalues[i])
        elseif type == "breath"
            Ht = total_Hamiltonian(H0,Hi,gvalues[i])
            H = total_Hamiltonian([Ht],Hb,omegas[i])
        end
        Energy[i] = real.(state' * H * state)
        dens[i,:] = density_profile_time(system,state,x=xs)
        norm[i] = sum(abs2.(state))
        fidelity[i] = abs.(state' * state_init)
    end
    results = (Energy,dens,norm,fidelity)
    return results
end


function density_profile_time(
    system :: Few_Particles_Hamonic_Oscillator,
    state;
    x :: Vector{Float64} = collect(-5.0:0.01:5.0)
    )

    particles = system.system.particles
    nho = system.system.lvls

    ndeg = size(particles, 1)
    OBDM = One_Body_Density_Matrix_time(system,state) 
    density_function = zeros(ComplexF64,(ndeg+2,length(x)))
    density_function[1,:] = x
    for p in 1:ndeg
        for i in 1:nho
            for j in 1:nho
                density_function[p+1, :] += OBDM[p, i, j] * (HO_wavefunctions(i-1, x) .* HO_wavefunctions(j-1, x))
            end
        end
    end
    density_function[end, :] = sum(density_function[2:end-1, :], dims=1)[1, :]
    return real.(density_function[end, :])
end

function One_Body_Density_Matrix_time(
    system :: Few_Particles_Hamonic_Oscillator,
    state
    )

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
    OBDM = zeros(ComplexF64,(ndeg,nho,nho))
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
    return OBDM
end