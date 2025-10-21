################################
#### Creates only the dictionary (to be implemented), loops on keys(dict)
#### It is possible to reduce the number of allocations?
################################

module BasisCreation_Develop
export do_basis

############## Optimized Code ##################
"""
    do_basis(particles, lvls, particletype) => basis, inbasis

Generate a basis for a quantum system based on the given particles, levels, and particle type.

# Inputs

    particles :: Vector{Int}            Vector specifying the number of particles in each internal state.
    lvls :: Int                         Number of energy levels available.
    particletype :: String              Type of particles, either `"fermi"` for fermions or `"bose"` for bosons.

# Returns

    basis :: Matrix{Int}                Matrix where each row represents a basis state.
    inbasis :: Dict{Vector{Int}, Int}   Dictionary mapping basis states to their indices.

# Notes
    The basis is generated with energy truncation [ref].

"""
function do_basis(
    particles :: Vector{Int},
    lvls :: Int,
    particletype :: String
    ) 

    N = sum(particles)

    if particletype == "fermi"
        lvlstates = lvls + div(sum(particles .* (particles .+ 1)), 2) - maximum(particles)
    elseif particletype == "bose"
        lvlstates = lvls
    else
        println("particletype: must be either 'fermi' or 'bose'")
        return
    end

    #basis = zeros(Int, numstates, N)
    globaliteration = 1
    #inbasis = Dict{Vector{Int}, Int}()
    basis = Dict{Vector{Int}, Int}()
    #basis, globaliteration, inbasis = iteration(basis, zeros(Int, N), particles, lvlstates, particletype, 1, globaliteration, inbasis)
    basis, globaliteration = iteration(basis, zeros(Int, N), particles, lvlstates, particletype, 1, globaliteration)
    return basis
    #return basis, inbasis
end

"""
    iteration_number_of_states_bose(i, j, particles, max, totalmax, states) => states

A helper function that calculates the number of states for a system of bosons with the energy truncation criteria.
"""
function iteration_number_of_states_bose(
    i :: Int, 
    j :: Int, 
    particles :: Vector{Int}, 
    max :: Int,
    totalmax :: Int, 
    states :: Int
    )

    for k in 1:max
        if i > 1
            states = iteration_number_of_states_bose(i-1,j,particles,min(k,totalmax-k+1),totalmax-k+1,states)
        elseif j > 1 
            states = iteration_number_of_states_bose(particles[j-1],j-1,particles,totalmax-k+1,totalmax-k+1,states)
        else
            states += 1
        end
    end
    return states
end

"""
    iteration_number_of_states_fermi(i, j, particles, max, totalmax, states) => states

A helper function that calculates the number of states for a system of fermions with the energy truncation criteria.
"""
function iteration_number_of_states_fermi(
    i :: Int, 
    j :: Int, 
    particles :: Vector{Int}, 
    max :: Int,
    totalmax :: Int, 
    states :: Int
    )
    
    while i > 0
        for k in i:max
            if i > 1
                ki = k + div(sum(particles[1:j-1] .* (particles[1:j-1] .+ 1)), 2) + div((i-1)*(i-2), 2) + 1
                states = iteration_number_of_states_fermi(i-1, j, particles, min(k-1, totalmax-ki+1), totalmax-k, states)
            elseif j > 1
                maxend = totalmax - div(sum(particles[1:j-2] .* (particles[1:j-2] .+ 1)), 2) - div(particles[j-1]*(particles[j-1]-1), 2) - k
                states = iteration_number_of_states_fermi(particles[j-1], j-1, particles, maxend, totalmax-k, states)
            else
                states += 1
            end
        end
        break
    end
    return states
end

""" 
    iteration() and creationparticle() are helper functions for generating the basis.

They are used to recursively create the basis states for the quantum system.

- 'iteration()' loops over the different particles (or internal states) and calls 'creationparticle()' 
    to generate the states for each particle type.
- 'creationparticle()' adds a particle in an allowed energy state and calls 'creationparticle()' again to add the next particle or 
        'iteration()' when there is no more particle of that type to add.
"""
function iteration(
    #basis :: Matrix{Int},
    basis :: Dict{Vector{Int}, Int},
    state :: Vector{Int},
    particles :: Vector{Int},
    lvls :: Int,
    particletype :: String,
    particlestate :: Int,
    globaliteration :: Int,
    #inbasis :: Dict{Vector{Int}, Int}
    ) 

    if particlestate <= length(particles)
        if particletype == "fermi"
            maxe = lvls
            maxnow = maxe - div(sum(particles[particlestate+1:end] .* (particles[particlestate+1:end] .+ 1)), 2) - div(particles[particlestate] * (particles[particlestate] - 1), 2)
            basis, globaliteration = creationparticle(basis, state, particles, maxe, particletype, particlestate, globaliteration, 0, maxnow)
        elseif particletype == "bose"
            basis, globaliteration = creationparticle(basis, state, particles, lvls, particletype, particlestate, globaliteration, 0, lvls)
        end
    else
        #basis[globaliteration, :] = state
        basis[state] = globaliteration
        globaliteration += 1
    end
    return basis, globaliteration
end

function creationparticle(
    #basis :: Matrix{Int},
    basis :: Dict{Vector{Int}, Int},
    state :: Vector{Int},
    particles :: Vector{Int}, 
    lvls :: Int,
    particletype :: String,
    particlestate :: Int,
    globaliteration :: Int,
    particleiteration :: Int,
    lastloc :: Int,
    #inbasis :: Dict{Vector{Int}, Int} 
    )

    if particles[particlestate] == 0
        return iteration(basis, state, particles, lvls, particletype, particlestate + 1, globaliteration)
    end
    if particletype == "fermi"
        for pos in (particles[particlestate] - particleiteration):lastloc
            statei = copy(state)
            statei[end] = (pos - 1) * length(particles) + particlestate
            statei = sort(statei,rev=true)
            #statei[(pos - 1) * length(particles) + particlestate] += 1
            if particleiteration + 1 != particles[particlestate]
                pose = pos + div(sum(particles[particlestate + 1:end] .* (particles[particlestate + 1:end] .+ 1)), 2) + div((particles[particlestate] - particleiteration - 1) * (particles[particlestate] - particleiteration - 2), 2)
                basis, globaliteration = creationparticle(
                    basis, statei, particles, lvls - pos, particletype, particlestate, globaliteration, particleiteration + 1, min(pos - 1, lvls - pose)
                )
            else
                basis, globaliteration, inbasis = iteration(
                    basis, statei, particles, lvls - pos, particletype, particlestate + 1, globaliteration
                )
            end
        end
    elseif particletype == "bose"
        for pos in 1:lastloc
            statei = copy(state)
            statei[end] = (pos - 1) * length(particles) + particlestate
            statei = sort(statei,rev=true)
            #statei[(pos - 1) * length(particles) + particlestate] += 1
            if particleiteration + 1 != particles[particlestate]
                basis, globaliteration = creationparticle(
                    basis, statei, particles, lvls - pos + 1, particletype, particlestate, globaliteration, particleiteration + 1, min(pos, lvls - pos + 1)
                )
            else
                basis, globaliteration = iteration(
                    basis, statei, particles, lvls - pos + 1, particletype, particlestate + 1, globaliteration
                )
            end
        end
    end

    return basis, globaliteration
end

end


module BasisCreationTotal_Develop
export do_basis

##############     Check!!     ##################
"""
    do_basis(particles, lvls, particletype) => basis, inbasis

Generate a basis for a quantum system based on the given particles, levels, and particle type.

# Inputs

    particles :: Vector{Int}            Vector specifying the number of particles in each internal state.
    lvls :: Int                         Number of energy levels available.
    particletype :: String              Type of particles, either `"fermi"` for fermions or `"bose"` for bosons.

# Returns

    basis :: Matrix{Int}                Matrix where each row represents a basis state.
    inbasis :: Dict{Vector{Int}, Int}   Dictionary mapping basis states to their indices.

# Notes
    The basis is generated without energy truncation [ref].

"""
function do_basis(
    particles :: Vector{Int},
    lvls :: Int,
    particletype :: String
    ) # no energy truncation.   

    internalstates = length(particles)
    if particletype=="fermi"
        numstates=1
        for i in 1:internalstates
            numstates=Int(numstates* binomial(lvls, particles[i]))
        end
        basis=zeros(Int,numstates,internalstates*lvls)
    elseif particletype=="bose"
        numstates=1
        for i in 1:internalstates
            numstates=Int(numstates*binomial(lvls+particles[i]-1,particles[i]))
        end
        basis=zeros(Int,numstates,internalstates*lvls)
    else
        println("Invalid particletype: must be 'fermi' or 'bose'")
        return
    end
    globaliteration=1
    inbasis=Dict{Vector{Int}, Int}() ## Use NTuple for better efficience?
    basis,globaliteration,inbasis=iteration(basis,zeros(Int,internalstates*lvls),particles,lvls,particletype,1,globaliteration,inbasis)
    return basis,inbasis
end

""" 
    iteration() and creationparticle() are helper functions for generating the basis.

They are used to recursively create the basis states for the quantum system.

- 'iteration()' loops over the different particles (or internal states) and calls 'creationparticle()' 
    to generate the states for each particle type.
- 'creationparticle()' adds a particle in an allowed energy state and calls 'creationparticle()' again to add the next particle or 
        'iteration()' when there is no more particle of that type to add.
"""
function iteration(
    basis :: Matrix{Int},
    state :: Vector{Int},
    particles :: Vector{Int},
    lvls :: Int,
    particletype :: String,
    particlestate :: Int,
    globaliteration :: Int,
    inbasis :: Dict{Vector{Int}, Int} 
    )

    if particlestate <= length(particles)
        basis,globaliteration=creationparticle(basis,state,particles,lvls,particletype,particlestate,globaliteration,0,lvls,inbasis)
    else
        basis[globaliteration,:]=state
        inbasis[state]=globaliteration
        globaliteration += 1
    end
    return basis,globaliteration,inbasis
end

function creationparticle(
    basis :: Matrix{Int},
    state :: Vector{Int},
    particles :: Vector{Int},
    lvls :: Int,
    particletype :: String,
    particlestate :: Int,
    globaliteration :: Int,
    particleiteration :: Int,
    lastloc :: Int,
    inbasis :: Dict{Vector{Int}, Int} 
    )
    
    if particles[particlestate]==0
        basis,globaliteration,inbasis=iteration(basis,state,particles,lvls,particletype,particlestate+1,globaliteration,inbasis)
    else
        if particletype=="fermi"
            for pos in (particles[particlestate]-particleiteration):lastloc
                statei=copy(state)
                statei[(pos-1)*length(particles)+particlestate] += 1
                if particleiteration+1 != particles[particlestate]
                    basis,globaliteration=creationparticle(basis,statei,particles,lvls,particletype,particlestate,globaliteration,particleiteration+1,pos-1,inbasis)
                else
                    basis,globaliteration,inbasis=iteration(basis,statei,particles,lvls,particletype,particlestate+1,globaliteration,inbasis)
                end
            end
        elseif particletype == "bose"
            for pos in 1:lastloc
                statei = copy(state)
                statei[(pos-1)*length(particles)+particlestate] += 1
                if particleiteration+1!=particles[particlestate]
                    basis,globaliteration=creationparticle(basis,statei,particles,lvls,particletype,particlestate,globaliteration,particleiteration+1,pos,inbasis)
                else
                    basis,globaliteration,inbasis=iteration(basis,statei,particles,lvls,particletype,particlestate+1,globaliteration,inbasis)
                end
            end
        end
    end
    return basis,globaliteration
end

end

