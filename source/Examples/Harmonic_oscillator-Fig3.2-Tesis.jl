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
particles = [1,1,1]
nho = 20
nstates = 50
type = "fermi"
#specific parameters for the calculation
target_state = 1
target_interaction = 20.0
glist = Vector(0:0.1:50)  
xlist = collect(-5.0:0.01:5.0)

println("Hamiltonian creation")
H_ho = HC.one_body_Hamiltonian(particles, nho, type)
H_ID = HC.two_body_Hamiltonian(particles, nho, type)

println("Diagonalization")
vals = Diag.energy_spectrum(H_ho, H_ID, glist, nstates, particles, nho)
vec, gef = Diag.eigenstate_calculation(H_ho, H_ID, target_interaction, target_state, particles, nho)

# the horizontal lines should be solved if we implement the feature of obtain the 
#    requested values of the interaction.
plot(title = "Energy spectrum",
     xlabel = "g",
     ylabel = "Energy",
     legend = nothing,
     xlims = (0, 20))
for i in 1:nstates
    plot!(vals[i,1,:], vals[i,2, :], label = "State $i", color=:blue)
end
display(plot!())

#Save the data
savefig("Energy_spectrum-Tesis3.2-particles3-nho20.png")
using JLD2 
jldsave("Energy_spectrum-Tesis3.2-particles3-nho20.jld2", true; vals, glist)



#Extra computations
"""
println("density calculation")
dens = Obs.density(vec, particles, nho, type,xlist)

println("OBDM calculation")
obdm_eig = Obs.One_Body_Density_Matrix_eigvals(vec, particles, nho, type)
obdm_space = Obs.One_Body_Density_Matrix_spatial(vec, particles, nho, type,xlist)

println("virial calculation")
virial = Obs.virial_energies(vec, particles, nho, type, gef)

#  Now, I store the xlist in the density variable. Mayve we can avoid that. 
#    If we use it as a structure (like dens.x , dens.rho, dens.rhototal ...)
plot(title = "density profile",
     xlabel = "x",
     ylabel = "Density",
     legend = :topright,
     xlims = (-5, 5),
     ylims = (0, 1))
for i in 1:2
    plot!(dens[1,:], dens[i+1, :], label = "State $i")
end
plot!(dens[1,:],dens[end,:],label="Total")
display(plot!())

# not so much to say here
plot(title = "One Body Density Matrix Eigenvalues",
     xlabel = "Eigenvalue index",
     ylabel = "Eigenvalue",
     legend = nothing,
     ylims = (0, 1))
scatter!(1:length(obdm_eig), obdm_eig, label = "Eigenvalues")
display(plot!())

#   Here I store the grid in the first index, and the sum on the last one. 
#     Maybe we can avoid that too.
plot(title = "One Body Density Matrix, particle 1",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,)
heatmap!(xlist, xlist, obdm_space[2,:,:], color=:blues, title="OBDM for state 1")
display(plot!())

plot(title = "One Body Density Matrix, particle 2",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,
     )
heatmap!(xlist, xlist, obdm_space[3,:,:], color=:blues, title="OBDM for state 2")
display(plot!())

"""
