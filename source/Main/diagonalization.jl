include("Sort.jl")

"""
    Computs the energy spectrum for the system.
    
        Input:  system -> The structure defines with Few_Particles_Hamonic_Oscillator, it have all the requiered information
                            It must to include the interaction values for whichyou want the spectrum in the correct format in "system.glist"
                            It must to include the number of states you want to compute in "system.numstates"
        Output: spectrum -> A structure that have two attributes, "spectrum.g" and "spectrum.energy"
                                "spectrum.g" is a three dimensional array, where the first index indicates the associated energy level (ground state, first excited...)
                                                                                 the second index is the array of interaction for that energy calculation
                                                                                 the third index is the value of the itneraction strength for each interaction pair
"""
function energy_spectrum(
    system :: Few_Particles_Hamonic_Oscillator
    )

    if length(system.H_mat.H0) == 0
        one_body_Hamiltonian!(system)
    end

    if length(system.H_mat.HI) == 0
        two_body_Hamiltonian!(system)
    end

    H0 = system.H_mat.H0
    Hi = system.H_mat.HI
    glist = system.glist
    num_states = system.numstates
    particles = system.system.particles
    sites = system.system.lvls

    sizeH = size(H0[1])[1]
    deg_g = size(Hi)[1]
    len_g = size(glist)[1]

    # If the format of glist is not the appropiate, here it tries to convert to the desired format
    if typeof(glist) == Matrix{Float64}
        if size(glist)[2] != deg_g
            println("Shape not correct. Computing SU(N) case with the first component values")
            temp = zeros((len_g,deg_g))
            for i in 1:deg_g
                temp[:,i] = glist[:,1]
            end
            glist = temp
        end
    else
        temp = zeros((len_g,deg_g))
        for i in 1:deg_g
            temp[:,i] = glist
        end
        glist = temp
    end
    totaleigenvals = zeros(len_g,num_states+deg_g)
    totaleigenvals[:,1:deg_g] = glist
    v0 = rand(sizeH)
    oldeigvecs = 1
    E0 = sum(particles)/2
    for (i,g) in enumerate(eachrow(glist))
        if ndims(g) == 0
            gt = [g]
        end
        gt = Vector(g)
        gcalc = zeros(size(gt)[1])

        # It calculates what should be th correction for the ground state using the previous energy calculation. 
        # With that, it stimates the value of the itneraction strenght that one should use to obtain the desired value after the correction.
        # This ensures to have the final result with almost the same range and intervals of the input grid of glist.
        for (k,gtem) in enumerate(gt)
            gcalc[k] = gtem/(1 + gtem * g_c(gtem,E0,sites,particles,only_gc_inv = true))
        end
        H = total_Hamiltonian(H0,Hi,gcalc)
        decomp, history = partialschur(
            H,
            v1 = v0,
            nev = num_states,
            which = :SR,
            tol = 1e-15,
        )
        # The system uses the supoerposition of the eigenstates of the prevous step to initialize the following calculation to reduce the computation time.
        v0 = vec(sum(decomp.Q; dims=2))
        eigvals = real.(decomp.eigenvalues)
        eigvecs = Matrix(decomp.Q)
        # Sorting protocol to have the eigenvalues sorted using an overlap criteria.
        if i != 1
            sort_eigs!(eigvecs,eigvals,oldeigvecs)
        end
        totaleigenvals[i,deg_g+1:end] = eigvals
        totaleigenvals[i,1:deg_g] = gcalc
        E0 = eigvals[1]
        
        oldeigvecs = deepcopy(eigvecs)

    end
    corrected_res = correction_eigenvalues(totaleigenvals,sites,particles,num_states)
    energies = totaleigenvals[:,deg_g+1:end]
    gvalues =  permutedims(corrected_res[:,1:deg_g,:],(1,3,2))
    spectrum = energy_spectrum_struct(gvalues,energies)
    return spectrum
end

"""
    Internal function to compute the eigenstate of a system for a given value of the interaction strength controlled by "system.target_g"
        It can be the ground state or some excitation, controlled by "system.target_state".
        The function computes the state corresponding to the physical value of the interaction strenght, it means that it is the
            state for that interaction after the correction. It could be time-consuming as i need to iterate to converge to the 
            correct numerical interaction strenght that returns the desired one after the correction.
            It also storr the numerical interaction strenght used in the calculation.
"""
function eigenstate!(
    system :: Few_Particles_Hamonic_Oscillator
    )

    if length(system.H_mat.H0) == 0
        one_body_Hamiltonian!(system)
    end

    if length(system.H_mat.HI) == 0
        two_body_Hamiltonian!(system)
    end

    H0 = system.H_mat.H0
    Hi = system.H_mat.HI
    g = system.target_g
    state = system.target_state
    particles = system.system.particles
    sites = system.system.lvls
    system.params = Internal_params((system.target_state, system.target_g),system.params.target_eigenstate,system.params.target_gcalc,system.params.OBDM,system.params.PCF)

    if ndims(g) == 0
        g = [g]
    end
    deg_g = size(Hi)[1]
    if size(g)[1] != deg_g
        g = g[1] .* ones(deg_g)
    end
    gef = g .+ 1
    gcal = copy(g)
    energy = 1.0
    sizeH = size(H0[1])[1]
    v0 = rand(sizeH)
    eigvecs = zeros(sizeH, state + 5)
    while abs(sum(gef .- g)) > 10^-10
        for k in 1:size(g)[1]
            gcal[k] = g[k] / (1+g[k]*g_c(g[k],energy,sites,particles,type = "fast",only_gc_inv = true))
        end
        H = total_Hamiltonian(H0,Hi,gcal)
        decomp, history = partialschur(
            H,
            v1 = v0,
            nev = state+5,
            which = :SR,
            tol = 1e-15,
        )
        v0 = vec(sum(decomp.Q; dims=2))
        eigvecs = Matrix(decomp.Q)
        eigvals = real.(decomp.eigenvalues)
        energy = eigvals[state]
        for k in 1:size(g)[1]
            gef[k] = g_c(gcal[k],energy,sites,particles,type = "fast")
        end
    end
    system.params = Internal_params(system.params.prev_params,eigvecs[:,state],system.params.target_gcalc,system.params.OBDM,system.params.PCF)
    system.params = Internal_params(system.params.prev_params,system.params.target_eigenstate,gcal,system.params.OBDM,system.params.PCF)
end

"""
    Internal function to correct all the interaction strenght values of the computed energy spectrum
"""
function correction_eigenvalues(
    vals :: Matrix{Float64},
    sites :: Int,
    particles :: Vector{Int},
    numstates :: Int; 
    type :: String = "fast"
    )

    size_vals = size(vals)
    deg_g = size_vals[2] - numstates
    numg = size_vals[1]
    corrected_results = zeros(numstates,deg_g+1,numg)
    for i in 1:numg
        for j in 1:numstates
            for k in 1:deg_g
                corrected_results[j,k,i] = g_c(vals[i,k],vals[i,j+deg_g],sites,particles,type = type)
            end
            corrected_results[j,deg_g+1,i] = vals[i,j+deg_g]
        end
    end
    return corrected_results
end

"""
    Internal function to compute the value of the correction for a given value of the energy and the interaction strength.
        It also depends on the system: number of particles and number of mode used.
"""
function g_c(
    g :: Float64,
    E :: Float64,
    sites :: Int,
    particles :: Vector{Int};
    type :: String = "fast",
    only_gc_inv :: Bool = false
    )

    N=sum(particles)
    nu = (E - N/2)/2
    M = floor( Int, (sites - 1) / 2 )
    if type == "fast"
        if nu == 0  
            nu = 1e-16
        end
        if nu < 0
            prefact = 1/(2 * pi * sqrt(2))
            first_part = 1/sqrt( -nu ) * atan( 2 * sqrt( (M + 1) * (-nu) ) / ( M + 1 + nu ) )
            gc_inv = (first_part) * prefact
        else
            prefact = 1/(2 * pi * sqrt(2))
            first_part = 1/sqrt(nu) * log( ( sqrt( M + 1 ) + sqrt( nu ) ) / ( sqrt( M + 1 ) - sqrt( nu ) ) )
            second_part = ( 1 + 1 /( 12 * ( M + 1 ) ) + 1 / ( 6 * ( M + 1 - nu ) ) ) / ( 2 * sqrt( M + 1 ) * ( M + 1 - nu) )
            gc_inv = (first_part + second_part) * prefact
        end
    elseif type == "precise" # seems that in some cases there are problems with the convergence
        prefact = 1/(2 * pi * sqrt(2))
        factor_1 = gamma( M + 3 / 2 ) * gamma( M - nu +1 )
        factor_2 = pFq( ( 1 , M + 3 / 2 , M - nu + 1 ) , ( M + 2 , M - nu + 2 ) , 1 ) / ( gamma( M + 2 ) * gamma( M - nu +2 ) )
        gc_inv = prefact * factor_1 * factor_2
    else
        println("type must be fast or precise")
        return g
    end
    if only_gc_inv
        return gc_inv
    end
    return g / ( 1 - g * gc_inv )
end