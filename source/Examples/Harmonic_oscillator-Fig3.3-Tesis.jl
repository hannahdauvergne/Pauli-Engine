### Example of use

####    The includes and imports is a temporal solution. 
####      All the usable functions may be in a single module?
include("../Develop/Hamiltonian_creation.jl")
include("../Develop/Diagonalization.jl")
include("../Develop/Observables.jl")

import .HamiltonianCreation_Develop as HC
import .Diagonalization_Develop as Diag
import .Observables_Develop as Obs

using Plots

# Parameters of the system. Coudl we create a structure for this?
particles = [1,1,1,1]
nho = 40
type = "fermi"
#specific parameters for the calculation
target_state = 1
target_interaction_vector = [0.0, 5.0, 10.0, 100.0]  # Vector of interactions to calculate
xlist = collect(-5.0:0.01:5.0)

println("Hamiltonian creation")
H_ho = HC.one_body_Hamiltonian(particles, nho, type)
H_ID = HC.two_body_Hamiltonian(particles, nho, type)


# Store the goundstates of the system for each target interaction
vec = Vector{Vector{Float64}}(undef, length(target_interaction_vector))
gef = Vector{Vector{Float64}}(undef, length(target_interaction_vector))
for (i, target_interaction) in enumerate(target_interaction_vector)
    println("Calculating for interaction g = $target_interaction")
    a, b = Diag.eigenstate_calculation(H_ho, H_ID, target_interaction, target_state, particles, nho)
    vec[i] = a 
    gef[i] = b
end

println("density calculation")
dens = Vector{Matrix{Float64}}(undef, length(target_interaction_vector))
for (i, target_interaction) in enumerate(target_interaction_vector)
    println("Calculating density for interaction g = $target_interaction")
    a = Obs.density(vec[i], particles, nho, type, xlist)
    dens[i] = a
end

#  Now, I store the xlist in the density variable. Mayve we can avoid that. 
#    If we use it as a structure (like dens.x , dens.rho, dens.rhototal ...)
plot(title = "density profile",
     xlabel = "x",
     ylabel = "Density",
     legend = :topright,
     xlims = (-3, 3),
     ylims = (0, 0.6))
for n in eachindex(target_interaction_vector)
    plot!(dens[n][1,:], dens[n][2, :], label="g = $(target_interaction_vector[n])")
end
#plot!(dens[1,:],dens[end,:],label="Total")
display(plot!())

#Save the data
savefig("GS_density-Tesis3.3-particles$(sum(particles))-nho40.png")
using JLD2 
jldsave("GS_density-Tesis3.3-particles$(sum(particles))-nho40.jld2", true; dens)

