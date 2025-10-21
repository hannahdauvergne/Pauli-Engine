##############
# Why seems more efficient in python the diagonalization as the basis increases? 
#HectorHere: can be because the algebra library beyond all functions is called with a different number of threads in python than in julia?
##############
module Diagonalization
### test which diagonalization method is the best...
###############   test with: nho = 20, particles = [1,1] (fermions), nstates = 11, glist = 0:0.1:10
###############     Arpak:                      343 ms      153764 allocations,     149.75 MiB
###############     Arpak Improved:             268 ms      121644 allocations,     148.85 MiB
###############     Arpak Parameters:           353 ms      139497 allocations,     154.36 MiB
###############     Arpak Parameters Improved:  250 ms       99011 allocations,     153.16 MiB       ## best for arpack. Implemented as the arpack option.
###############     Arnoldi:                    265 ms       35643 allocations,     149.67 MiB
###############     Arnoldi Improved:           263 ms       34734 allocations,     141.47 MiB
###############     Arnoldi Improved v2:        226 ms       32078 allocations,     149.34 MiB       ## best for arnoldi (and in general. Implemented as the arnoldi option, and the default one.
###############     Scipy (Python):             782 ms      
###############
###############   test with: nho = 40, particles = [1,1] (fermions), nstates = 11, glist = 0:0.1:10
###############     Arpak:                      6.076 s       245180 allocations,    2.05 GiB
###############     Arpak Improved:             4.777 s       190799 allocations,    2.05 GiB
###############     Arpak Parameters:           5.626 s       214417 allocations,    2.07 GiB
###############     Arpak Parameters Improved:  4.499 s       167678 allocations,    2.06 GiB
###############     Arnoldi:                    4.501 s        50426 allocations,    2.05 GiB  
###############     Arnoldi Improved:           4.536 s        49791 allocations,    2.02 GiB
###############     Arnoldi Improved v2:        4.028 s        45141 allocations,    2.05 GiB
###############     Scipy (Python):             8.168 s 
###############
###############   test with: nho = 60, particles = [1,1] (fermions), nstates = 11, glist = 0:0.1:10
###############     Arpak:                      60.024 s    314884 allocations,    10.09 GiB
###############     Arpak Improved:             48.550 s    251967 allocations,    10.09 GiB
###############     Arpak Parameters:           55.988 s    266033 allocations,    10.12 GiB
###############     Arpak Parameters Improved:  43.617 s    217774 allocations,    10.12 GiB
###############     Arnoldi:                    45.877 s     62359 allocations,    10.10 GiB  
###############     Arnoldi Improved:           45.953 s     61718 allocations,    10.03 GiB
###############     Arnoldi Improved v2:        40.183 s     56071 allocations,    10.10 GiB
###############     Scipy (Python):             55.878 s
using Arpack
using KrylovKit # not implemented yet
using ArnoldiMethod
using SparseArrays
using LinearAlgebra # I dont know if it is used
using HypergeometricFunctions # used in the correction functions, but it is so slow
using SpecialFunctions          # and have large numerical errors. maybe we should not use it... 
include("Sort.jl")

export energy_spectrum,eigenstate_calculation

##### Sort works better if g starts at 0
######## Maybe implementation of computing first repulsive from 0 and then atractive from 0
##  It need particles and nho, maybe we can avoid that
"""
    energy_spectrum(H0, Hi, glist, num_states, particles, sites; method, sorted, corrected) => totaleigenvals
Computes the lower eigenvalues of a system for several interaction strengths.

# Inputs

    H0 :: Vector{SparseMatrixCSC{Float64, Int}}     Single-particle Hamiltoinan
    Hi :: Vector{SparseMatrixCSC{Float64, Int}}     Interaction Hamiltoinan
    glist :: Vector{Float64}                        Values of the interaction to compute the eigenvalues
    num_states :: Int                               Number of eigenvalues requiered
    particles :: Vector{Int}                        Number of particles in the system    
    sites :: Int                                    Number of sites used in the calculation

# Optional parameters
    method :: String = "arnoldi"                    Eigensolver method, "arnoldi" or "arpack"
    sorted :: Bool = true                           Use of a routine to have the eigenstates sorted by overlap criteria
    corrected :: Bool = true                        Use of a correction on the calculation due to the finite number of modes [ref]

# Returns
    
    totaleigenvals :: Matrix{Float64,Float64}       Depends on if it is corrected or not. To be improved
"""
function energy_spectrum(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    glist :: Union{Vector{Float64}, Matrix{Float64}},
    num_states :: Int,
    particles :: Vector{Int},
    sites :: Int;
    method :: String = "arnoldi",
    sorted :: Bool = true,
    corrected :: Bool = true,
    request :: Bool = false,
    fast :: Bool = true
    )

    if request
        if fast
            totaleigenvals = eig_arnoldi_improved_v2_sorted_corrected_fast(H0, Hi, glist, num_states,particles,sites)
        else
            totaleigenvals = eig_arnoldi_improved_v2_sorted_corrected_complete(H0, Hi, glist, num_states,particles,sites)
        end
    else

        if method == "arnoldi"
            if sorted
                totaleigenvals = eig_arnoldi_improved_v2_sorted(H0, Hi, glist, num_states)
            else
                totaleigenvals = eig_arnoldi_improved_v2(H0, Hi, glist, num_states)
            end
        elseif method == "arpack"
            if sorted
                totaleigenvals = eigenvalues_params_improved_sorted(H0, Hi, glist, num_states)
            else
                totaleigenvals = eigenvalues_params_improved(H0, Hi, glist, num_states)
            end
        else
            error("Method must be 'arnoldi' or 'arpack'")
        end
        if corrected
            totaleigenvals = correction_eigenvalues(totaleigenvals, sites, particles,num_states, type="fast")
        end
    end
    return totaleigenvals
end

"""
    eigensstate_calculation(H0, Hi, g, state, particles, sites; method, corrected) => eigenstate, gcalc
Computes a given eigenstate for a given value of the itneraction.

# Inputs 

    H0 :: Vector{SparseMatrixCSC{Float64, Int}}     Single-particle Hamiltoinan
    Hi :: Vector{SparseMatrixCSC{Float64, Int}}     Interaction Hamiltoinan
    g :: Float64                                    Value of the interaction to compute the eigenstate
    state :: Int                                    Index of the eigenstate requeired
    particles :: Vector{Int}                        Number of particles in the system    
    sites :: Int                                    Number of sites used in the calculation

# Optional parameters
    method :: String = "arnoldi"                    Eigensolver method, "arnoldi" or "arpack"
    corrected :: Bool = true                        Use of a correction on the calculation due to the finite number of modes [ref]

# Returns

    eigenstate :: Vector{Float64}                   Eigenstate requiered for the physical interaction strength on the input.
    gcalc :: Vector{Float64}                        Value of the interaction used in the numerical calulation.

"""
function eigenstate_calculation(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    g :: Float64,
    state :: Int,
    particles :: Vector{Int},
    sites :: Int;
    method :: String = "arnoldi",
    corrected :: Bool = true
    )
    if method == "arnoldi"
        if corrected
            eigenstate,gcalc = state_arnoldi_correction(H0, Hi, g, state, particles, sites)
        else
            eigenstate = state_arnoldi(H0, Hi, g, state)
            gcalc = g .* ones(size(Hi, 1))
        end
    elseif method == "arpack"
        if corrected
            eigenstate,gcalc = state_arpack_correction(H0, Hi, g, state, particles, sites)
        else
            eigenstate = state_arpack(H0, Hi, g, state)
            gcalc = g .* ones(size(Hi, 1))
        end
    else
        error("Method must be 'arnoldi' or 'arpack'")
    end
    return eigenstate,gcalc
end

function eig_arnoldi(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    glist :: Vector{Float64},
    num_states :: Int
    )

    len_g = size(glist)[1]
    totaleigenvals = zeros(len_g,num_states+1)
    totaleigenvals[:,1] = glist
    for (i,g) in enumerate(glist)
        if ndims(g) == 0
            g = [g]
        end
        H = total_Hamiltonian(H0,Hi,g)
        decomp, history = partialschur(
            H,
            nev = num_states,
            which = :SR,
            tol = 1e-15
        )
        totaleigenvals[i,2:end] = real.(decomp.eigenvalues)
    end
    return totaleigenvals
end


### In principle it should converge faster, as it uses information from the previous step
### but seems that needs the same amount of iterations
function eig_arnoldi_improved(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    glist :: Vector{Float64},
    num_states :: Int
    )

    sizeH = size(H0[1])[1]
    len_g = size(glist)[1]
    totaleigenvals = zeros(len_g,num_states+1)
    totaleigenvals[:,1] = glist
    V_in, H_in = rand(sizeH, max(21,1+2*num_states)), rand(max(1+2*num_states,21), max(2*num_states,20))
    arnoldi = ArnoldiWorkspace(V_in, H_in)
    for (i,g) in enumerate(glist)
        if ndims(g) == 0
            g = [g]
        end
        H = total_Hamiltonian(H0,Hi,g)
        decomp, history = partialschur!(
            H,
            arnoldi,
            nev = num_states,
            which = :SR,
            tol = 1e-15,
        )
        totaleigenvals[i,2:end] = real.(decomp.eigenvalues)
    end
    return totaleigenvals
end

function eig_arnoldi_improved_v2(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    glist :: Vector{Float64},
    num_states :: Int
    )

    sizeH = size(H0[1])[1]
    len_g = size(glist)[1]
    totaleigenvals = zeros(len_g,num_states+1)
    totaleigenvals[:,1] = glist
    v0 = rand(sizeH)
    for (i,g) in enumerate(glist)
        if ndims(g) == 0
            g = [g]
        end
        H = total_Hamiltonian(H0,Hi,g)
        decomp, history = partialschur(
            H,
            v1 = v0,
            nev = num_states,
            which = :SR,
            tol = 1e-15,
        )
        v0 = vec(sum(decomp.Q; dims=2))
        totaleigenvals[i,2:end] = real.(decomp.eigenvalues)
    end
    return totaleigenvals
end

function eig_arnoldi_improved_v2_sorted(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    glist :: Union{Vector{Float64},Matrix{Float64}},
    num_states :: Int
    )

    sizeH = size(H0[1])[1]
    len_g = size(glist)[1]
    if typeof(glist) == Matrix{Float64}
        deg_g = size(glist)[2]
    else
        deg_g = 1
    end
    totaleigenvals = zeros(len_g,num_states+deg_g)
    totaleigenvals[:,1:deg_g] = glist
    v0 = rand(sizeH)
    oldeigvecs = 1
    for (i,g) in enumerate(eachrow(glist))
        if ndims(g) == 0
            g = [g]
        end
        g = Vector(g)
        H = total_Hamiltonian(H0,Hi,g)
        decomp, history = partialschur(
            H,
            v1 = v0,
            nev = num_states,
            which = :SR,
            tol = 1e-15,
        )
        v0 = vec(sum(decomp.Q; dims=2))
        eigvals = real.(decomp.eigenvalues)
        eigvecs = Matrix(decomp.Q)
        if i != 1
            sort_eigs!(eigvecs,eigvals,oldeigvecs)
        end
        totaleigenvals[i,deg_g+1:end] = eigvals
        
        oldeigvecs = deepcopy(eigvecs)

    end
    return totaleigenvals
end

function eig_arnoldi_improved_v2_sorted_corrected_fast(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    glist :: Union{Vector{Float64},Matrix{Float64}},
    num_states :: Int,
    particles :: Vector{Int},
    sites :: Int
    )

    sizeH = size(H0[1])[1]
    len_g = size(glist)[1]
    if typeof(glist) == Matrix{Float64}
        deg_g = size(glist)[2]
    else
        deg_g = 1
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
        v0 = vec(sum(decomp.Q; dims=2))
        eigvals = real.(decomp.eigenvalues)
        eigvecs = Matrix(decomp.Q)
        if i != 1
            sort_eigs!(eigvecs,eigvals,oldeigvecs)
        end
        totaleigenvals[i,deg_g+1:end] = eigvals
        E0 = eigvals[1]
        
        oldeigvecs = deepcopy(eigvecs)

    end
    return totaleigenvals
end


function eig_arnoldi_improved_v2_sorted_corrected_complete(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    glist :: Union{Vector{Float64},Matrix{Float64}},
    num_states :: Int,
    particles :: Vector{Int},
    sites :: Int
    )

    sizeH = size(H0[1])[1]
    len_g = size(glist)[1]
    if typeof(glist) == Matrix{Float64}
        deg_g = size(glist)[2]
    else
        deg_g = 1
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
        for (k,gtem) in enumerate(gt)
            gcalc[k] = gtem/(1 + gtem * g_c(gtem,E0,sites,particles,only_gc_inv = true))
        end
        convergedvals = zeros(num_states)
        eigvals = zeros(num_states)
        convergedvecs = zeros(sizeH,num_states)
        eigvecs = zeros(sizeH,num_states)
        for j in 1:num_states
            gprev = gcalc .+ 1
            while abs(sum(gprev .- gt)) > 1e-10
                H = total_Hamiltonian(H0,Hi,gcalc)
                decomp, history = partialschur(
                    H,
                    v1 = v0,
                    nev = num_states,
                    which = :SR,
                    tol = 1e-15,
                )
                v0 = vec(sum(decomp.Q; dims=2))
                eigvals = real.(decomp.eigenvalues)
                for (nn,gp) in enumerate(gcalc)
                    gprev[nn] = gp / (1 - gp * g_c(gp,eigvals[j],sites,particles,only_gc_inv = true))
                end
                eigvecs = Matrix(decomp.Q)
                for (k,gtem) in enumerate(gt)
                    gcalc[k] = gtem/(1 + gtem * g_c(gtem,eigvals[j],sites,particles,only_gc_inv = true))
                end
            end
            convergedvals[j] = copy(eigvals[j])
            convergedvecs[:,j] = copy(eigvecs[:,j])
            E0 = sum(convergedvals) / num_states
        end
        if i != 1
            sort_eigs!(convergedvecs,convergedvals,oldeigvecs)
        end
        totaleigenvals[i,deg_g+1:end] = convergedvals
        E0 = eigvals[1]
        
        oldeigvecs = deepcopy(convergedvecs)

    end
    return totaleigenvals
end

function state_arnoldi(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    g :: Float64,
    state :: Int
    )

    if ndims(g) == 0
        g = [g]
    end
    H = total_Hamiltonian(H0,Hi,g)
    decomp, history = partialschur(
        H,
        nev = state + 5,
        which = :SR,
        tol = 1e-15
    )
    return Matrix(decomp.Q)[:,state]
end

function state_arnoldi_correction(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    g :: Float64,
    state :: Int,
    particles :: Vector{Int},
    sites :: Int;
    tol :: Float64 = 1e-10
    )

    if ndims(g) == 0
        g = [g]
    end
    gef = g .+ 1
    gcal = copy(g)
    energy = 1.0
    sizeH = size(H0[1])[1]
    v0 = rand(sizeH)
    eigvecs = zeros(sizeH, state + 5)
    while abs(sum(gef .- g)) > tol
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
    return eigvecs[:,state], gcal
end

function eigenvalues(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    glist :: Vector{Float64},
    num_states :: Int
    )

    len_g = size(glist)[1]
    totaleigenvals = zeros(len_g,num_states+1)
    totaleigenvals[:,1] = glist
    for (i,g) in enumerate(glist)
        if ndims(g) == 0
            g = [g]
        end
        H = total_Hamiltonian(H0,Hi,g)
        vals,vecs = eigs(H,
            nev = num_states,
            which = :SR,
            maxiter = 5000,
            tol = 0,
            ncv =  2*num_states+1
        )
        totaleigenvals[i,2:end] = vals
    end
    return totaleigenvals
end

function eigenvalues_params(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    glist :: Vector{Float64},
    num_states :: Int
    )

    D = size(H0[1],1)
    len_g = size(glist)[1]
    totaleigenvals = zeros(len_g,num_states+1)
    totaleigenvals[:,1] = glist
    for (i,g) in enumerate(glist)
        if ndims(g) == 0
            g = [g]
        end
        H = total_Hamiltonian(H0,Hi,g)
        vals,vecs = eigs(H,
            nev = num_states,
            which = :SR,
            maxiter = 10*D,
            tol = 0,
            ncv = 2*min(D,max(20,2*num_states+1))
        )
        totaleigenvals[i,2:end] = vals
    end
    return totaleigenvals
end

function eigenvalues_improved(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    glist :: Vector{Float64},
    num_states :: Int
    )

    sizeH = size(H0[1])[1]
    len_g = size(glist)[1]
    totaleigenvals = zeros(len_g,num_states+1)
    totaleigenvals[:,1] = glist
    v0_t = rand(sizeH)
    for (i,g) in enumerate(glist)
        if ndims(g) == 0
            g = [g]
        end
        H = total_Hamiltonian(H0,Hi,g)
        vals,vecs = eigs(H,
            nev = num_states,
            which = :SR,
            maxiter = 5000,
            tol = 0,
            ncv =  2*num_states+1,
            v0 = v0_t
        )
        totaleigenvals[i,2:end] = vals
        v0_t = sum(vecs; dims=2)
        v0_t = vec(v0_t) 
    end
    return totaleigenvals
end

function eigenvalues_params_improved(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    glist :: Vector{Float64},
    num_states :: Int
    )

    D = size(H0[1],1)
    len_g = size(glist)[1]
    totaleigenvals = zeros(len_g,num_states+1)
    totaleigenvals[:,1] = glist
    v0_t = rand(D)
    for (i,g) in enumerate(glist)
        if ndims(g) == 0
            g = [g]
        end
        H = total_Hamiltonian(H0,Hi,g)
        vals,vecs = eigs(H,
            nev = num_states,
            which = :SR,
            maxiter = 10*D,
            tol = 0,
            ncv = 2*min(D,max(20,2*num_states+1)),
            v0 = v0_t
        )
        totaleigenvals[i,2:end] = vals
        v0_t = vec(sum(vecs; dims=2))
    end
    return totaleigenvals
end

function eigenvalues_params_improved_sorted(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    glist :: Vector{Float64},
    num_states :: Int
    )

    D = size(H0[1],1)
    len_g = size(glist)[1]
    totaleigenvals = zeros(len_g,num_states+1)
    totaleigenvals[:,1] = glist
    v0_t = rand(D)
    oldeigvecs = 1
    for (i,g) in enumerate(glist)
        if ndims(g) == 0
            g = [g]
        end
        H = total_Hamiltonian(H0,Hi,g)
        vals,vecs = eigs(H,
            nev = num_states,
            which = :SR,
            maxiter = 10*D,
            tol = 0,
            ncv = 2*min(D,max(20,2*num_states+1)),
            v0 = v0_t
        )
        if i != 1
            sort_eigs!(vecs,vals,oldeigvecs)
        end
        totaleigenvals[i,2:end] = vals
        oldeigvecs = deepcopy(vecs)
        v0_t = vec(sum(vecs; dims=2))
    end
    return totaleigenvals
end

function state_arpack(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    g :: Float64,
    state :: Int
    )

    if ndims(g) == 0
        g = [g]
    end
    H = total_Hamiltonian(H0,Hi,g)
    vals,vecs = eigs(H,
        nev = state + 5,
        which = :SR,
        maxiter = 5000,
        tol = 0,
        ncv = 2*state+1
    )
    return vecs[:,state]
end

function state_arpack_correction(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    g :: Float64,
    state :: Int,
    particles :: Vector{Int},
    sites :: Int;
    tol :: Float64 = 1e-10
    )

    if ndims(g) == 0
        g = [g]
    end
    gef = g .+ 1
    gcal = copy(g)
    energy = 1
    D = size(H0[1],1)
    v0_t = rand(D)
    while abs(sum(gef .- g)) > tol
        for k in 1:size(g)[1]
            gcal[k] = g[k] / (1+g[k]*g_c(g[k],energy,sites,particles,type = "fast",only_gc_inv = true))
        end
        H = total_Hamiltonian(H0,Hi,gcal)
        vals,vecs = eigs(H,
            nev = state + 5,
            which = :SR,
            maxiter = 10*D,
            tol = 0,
            ncv = 2*min(D,max(20,2*(state+5)+1)),
            v0 = v0_t
        )
        v0_t = vec(sum(vecs; dims=2))
        eigvecs = Matrix(vecs)
        eigvals = real.(vals)
        energy = eigvals[state]
        for k in 1:size(g)[1]
            gef[k] = g_c(gcal[k],energy,sites,particles,type = "fast")
        end
    end

    return eigvecs[:,state],gcal
end

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

function total_Hamiltonian(
    H0 :: Vector{SparseMatrixCSC{Float64, Int}},
    Hi :: Vector{SparseMatrixCSC{Float64, Int}},
    g :: Vector{Float64}
    )

    ideg_int = size(Hi)

    if length(g) != ideg_int[1]
        g = g[1] * ones(ideg_int[1])
    end

    H= sum(H0)+sum(g .*Hi)
    return H
end

end