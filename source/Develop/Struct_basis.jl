"""
    Creates the many body basis for a given system, by considering a truncation using an energy cuttoff [ref]
"""
function basis_creation(
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

    globaliteration = 1
    basis = Dict{Vector{Int}, Int}()
    basis, globaliteration = iteration(basis, zeros(Int, N), particles, lvlstates, particletype, 1, globaliteration)
    return basis
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
    basis :: Dict{Vector{Int}, Int},
    state :: Vector{Int},
    particles :: Vector{Int},
    lvls :: Int,
    particletype :: String,
    particlestate :: Int,
    globaliteration :: Int,
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
        basis[state] = globaliteration
        globaliteration += 1
    end
    return basis, globaliteration
end

function creationparticle(
    basis :: Dict{Vector{Int}, Int},
    state :: Vector{Int},
    particles :: Vector{Int}, 
    lvls :: Int,
    particletype :: String,
    particlestate :: Int,
    globaliteration :: Int,
    particleiteration :: Int,
    lastloc :: Int,
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
                basis, globaliteration = iteration(
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

function basis_parity(    
    particles :: Vector{Int},
    lvls :: Int,
    particletype :: String,
    parity :: Bool
    )

    basis_t = basis_creation(particles, lvls, particletype)

    deg = length(particles)

    basis = Dict()
    k_indx = 1
    for key in keys(basis_t)
        state_parity = sum(div.(key.-1 ,deg)) % 2
        if state_parity != parity 
            basis[key] = k_indx
            k_indx += 1
        end
    end

    return basis
end