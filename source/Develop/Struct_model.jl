"""
    Main information of the system. In order to avoid problems, it should not changed after the initialization.
    It includes the number of particles, the number of modes used and the type of the particles (bosons or fermions)
    Also includes the many-body basis (so, the system only need to computes once)
    and an id code to identify the system, that can be useful if one want to save the files.
"""
struct System
    particles :: Vector{Int}
    lvls :: Int
    particletype :: String
    basis :: Dict{Vector{Int}, Int}
    id :: String
end

"""
    Structure to save the Hamiltonian matrices, it should only be modified by the creation hamiltonian functions to 
    not store Hamiltonians that does not correspond to the system.
"""
struct Hamiltonians
    H0 :: Vector{SparseMatrixCSC{Float64, Int64}}
    HI :: Vector{SparseMatrixCSC{Float64, Int64}}
    Hkin :: Vector{SparseMatrixCSC{Float64, Int64}}
    Hpot :: Vector{SparseMatrixCSC{Float64, Int64}}
end

"""
    Internal parameters of the system that does not should be modified except for the corresponding functions
        Stores information as the analyzed eigenstate, its numerical interaction strenght, the OBDM.
            as are things that could be needed more than one time and with that it can be the code more efficient.
        Also stores the previous parameters of the calculations, if there is the need to compute these things again for a different parameters.
"""
struct Internal_params
    prev_params :: Tuple{Int,Vector{Float64} }
    target_eigenstate :: Vector{Float64}
    target_gcalc :: Vector{Float64}
    OBDM :: Array{Float64, 3}
    PCF :: Array{Float64, 6} 
end

"""
    Main structure. It includes the system information, the Hamiltonians and the internal parameters. 
    It also have the editable attributes to compute the desired observables. 
    The values of ht einteraction strength for the energy spectrum as well the number of states computed.
    The state and the itneraction on what you want to compute the density, OBM and so.
    The path where one wants to save the files, if want.
"""
mutable struct Few_Particles_Hamonic_Oscillator  
    system :: System
    H_mat :: Hamiltonians
    params :: Internal_params
    glist :: Union{Vector{Float64},Matrix{Float64}}
    numstates :: Int
    target_state :: Int
    target_g :: Union{Int,Float64,Vector{Float64}}
    path :: String
    path_H :: String
    save_local :: Bool
end

"""
    Structure to store the density profile information.
    It contains the spatial grid in "x"
    The densiti of each component in "dens"
    the total density in "dentotal"
    and the value of the itneraction strenght used in the numerical calulation in "g". This does not correspond to the physical one.
"""
struct density_struct
    x :: Vector{Float64}
    dens :: Matrix{Float64}
    denstotal :: Vector{Float64}
    g :: Vector{Float64}
end

"""
    Structure to store the results of the energy spectrum
    Contains the final result of the physical valu of interaction strenght for each eigenstate an all the computed points in "g"
    Contains the energy for each eigenstate and interaction computed in "energy"
"""
struct energy_spectrum_struct
    g :: Array{Float64, 3}
    energy :: Matrix{Float64}
end

"""
    Inizialitacion function. As inputs you must provide the number of particles.
        If you not introduce the number of h.o. levels, it asing to a default value (not recomended)
        If you not specify the nature of the particles, it assumes that are bosons.
    Apart of the basic system information, you can initilize the additional parametes of the object.
    If are not initialized, it uses a deaul values.
    These elements can be modified after the inizilitacion.
"""
function Few_Particles_Hamonic_Oscillator(
    particles,
    lvls=nothing,
    particletype="bose";
    glist = collect(0:0.1:15),
    numstates = 2,
    target_state = 1,
    target_g = 0,
    path = "",
    path_H = "Hamiltonians/",
    save_local = false,
    id = "default",
    flip_component = false,
    parity_restriction = false,
    parity = true
    )

    H0 = Vector{SparseMatrixCSC{Float64, Int64}}()
    HI = Vector{SparseMatrixCSC{Float64, Int64}}()
    Hkin = Vector{SparseMatrixCSC{Float64, Int64}}()
    Hpot = Vector{SparseMatrixCSC{Float64, Int64}}()
    prev_params = (0,Vector{Float64}())
    target_eigenstate = Vector{Float64}(undef,0)
    target_gcalc = Vector{Float64}()
    OBDM = Array{Float64,3}(undef, 0, 0, 0)
    PCF = Array{Float64,6}(undef, 0, 0, 0, 0, 0, 0)

    if isnothing(lvls)
        lvls = Int(min(100,ceil(50000^(1/sum(particles))))) # Find a good authomatic value. This formula is temporal
        println("using $lvls states")
    end

    if flip_component
        if id == "default"
            id = particletype * "_" * string(sum(particles)) * "_SOC_nho_" * string(lvls)
        end
        particles = sum(particles)

        basis = Dict{Vector{Int}, Int}()
        for part in 0:particles
            temp_particles = particles - part
            temp_basis = basis_creation([temp_particles, part], lvls, particletype)
            offset = length(basis)
            temp_dic = Dict(k => v + offset for (k, v) in temp_basis)
            basis = merge(basis, temp_dic)
        end
        particles = [particles,0]
    elseif parity_restriction
        if id == "default"
            if parity
                id = particletype * "_" * join(string.(particles), "_") *"_nho_" * string(lvls) * "even_parity"
            else
                id = particletype * "_" * join(string.(particles), "_") *"_nho_" * string(lvls) * "odd_parity"
            end
        end
        basis = basis_parity(particles, lvls,particletype,parity)
    else
        if id == "default"
            id = particletype * "_" * join(string.(particles), "_") * "_nho_" * string(lvls)
        end

        basis = basis_creation(particles, lvls,particletype)
    end


    new_system = Few_Particles_Hamonic_Oscillator(
        System(particles,lvls,particletype,basis,id),
        Hamiltonians(H0,HI,Hkin,Hpot),
        Internal_params(prev_params,target_eigenstate,target_gcalc,OBDM,PCF),
        glist,numstates,
        target_state,target_g,
        path,path_H,save_local
        )
    return new_system
end

