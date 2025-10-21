### Example of use

####    The includes and imports is a temporal solution. 
####      All the usable functions may be in a single module?
include("../Develop/Basis_creation.jl")
include("../Develop/Hamiltonian_creation.jl")
include("../Develop/Diagonalization.jl")
include("../Develop/Observables.jl")

import .BasisCreation_Develop as BC
import .HamiltonianCreation_Develop as HC
import .Diagonalization_Develop as Diag
import .Observables_Develop as Obs

using Plots

# Parameters of the system. Coudl we create a structure for this?
particles = [1,1,1,1]
nho = 20
type = "fermi"
#specific parameters for the calculation
interaction_strength = 55.0
target_state = 1
xlist = collect(-5.0:0.01:5.0)

println("Hamiltonian creation")
H_ho = HC.one_body_Hamiltonian(particles, nho, type)
H_ID = HC.two_body_Hamiltonian(particles, nho, type)

"""
# Order of interaction strengths [1-1, 1-2, 1-3, ... 1-n, 2-2, 2-3, ... 2-n, 3-3, ... n-n]
# As they are fermions, the j-j interaction will not matter, so we can set it to zero.
"""


"""
First system: Connectivity interaction graph
A --- B 
|     |    
|     |
C     D
"""
target_interaction = interaction_strength * [0, 1, 1, 0, 0, 1, 0, 0, 0, 0] 

H = Diag.total_Hamiltonian(H_ho, H_ID, target_interaction)

# Here it breaks because the eigenstate calculation function does not accept the interaction vector, only a float
#TODO: Fix this to allow many interactions at once
# Store the goundstates of the system 
fec, gef = Diag.eigenstate_calculation(H_ho, H_ID, target_interaction, target_state, particles, nho)

println("density calculation")
dens = Obs.density(vec[i], particles, nho, type, xlist)

#  Now, I store the xlist in the density variable. Maybe we can avoid that. 
#    If we use it as a structure (like dens.x , dens.rho, dens.rhototal ...)
plot(title = "density profile",
     xlabel = "x",
     ylabel = "Density",
     legend = :topright,
     xlims = (-4, 4),
     ylims = (0, 0.8))
for n in eachindex(particles)
    plot!(dens[1,:], dens[n+1, :], label="Specie $n")
end
#plot!(dens[1,:],dens[end,:],label="Total")
display(plot!())

#Save the data
savefig("GS_density-Tesis3.5-g$target_interaction-nho$nho.png")
using JLD2 
jldsave("GS_density-Tesis3.5-g$target_interaction-nho$nho.jld2", true; dens)
