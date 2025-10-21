### Example of use

####    The includes and imports is a temporal solution. 
####      All the usable functions may be in a single module?
include("Backup/Basis_creation.jl")
include("Backup/Hamiltonian_creation.jl")
include("Backup/Diagonalization.jl")
include("Backup/Observables.jl")

import .BasisCreation as BC
import .HamiltonianCreation as HC
import .Diagonalization as Diag
import .Observables as Obs

using Plots

# Parameters of the system. Coudl we create a structure for this?
particles = [1, 1]
nho = 10
nstates = 5
type = "bose"
#specific parameters for the calculation
target_state = 1
target_interaction = 20.0
glist = Vector(0:0.1:50)  # this is the range of g prior correction. We should adapt the calculation to return this range after correction?
xlist = collect(-5.0:0.01:5.0)

println("basis creation")
basis, dic_base = BC.do_basis(particles,nho,type)

println("Hamiltonian creation")
H_ho = HC.one_body_Hamiltonian(particles, nho, type)
H_ID = HC.two_body_Hamiltonian(particles, nho, type)

println("Diagonalization")
vals = Diag.energy_spectrum(H_ho, H_ID, glist, nstates, particles, nho)
vec, gef = Diag.eigenstate_calculation(H_ho, H_ID, target_interaction, target_state, particles, nho)

println("density calculation")
dens = Obs.density(vec, particles, nho, type,xlist)

println("OBDM calculation")
obdm_eig = Obs.One_Body_Density_Matrix_eigvals(vec, particles, nho, type)
obdm_space = Obs.One_Body_Density_Matrix_spatial(vec, particles, nho, type,xlist)

println("virial calculation")
virial = Obs.virial_energies(vec, particles, nho, type, gef)

# the horizontal lines should be solved if we implement the feature of obtain the 
#    requested values of the itneraction.
plot(title = "Energy spectrum",
     xlabel = "g",
     ylabel = "Energy",
     legend = nothing,
     xlims = (0, 50))
for i in 1:nstates
    plot!(vals[i,1,:], vals[i,2, :], label = "State $i")
end
display(plot!())

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

