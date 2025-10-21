include("Main/main.jl")
using .Harmonic_Oscillator_solver
using Plots
using LaTeXStrings

# set the number of particles as a mandatory input
particles = [1,1,1] 
# creates the system struct. It is possible to modify any parameter after that, but be careful as after a change some parameters would not be updated. Eg, if you change the number of particles, the basis won't be authomatically updated.
system = Few_Particles_Hamonic_Oscillator(particles,10, "bose")
system2 = Few_Particles_Hamonic_Oscillator(particles)
system3 = Few_Particles_Hamonic_Oscillator([2,3],15,"fermi")
system4 = Few_Particles_Hamonic_Oscillator([2,2],15)

# you can compute directly observables. If some parameter must be computed (Hamiltonian for example) It does only once
# Maybe we can add a way to store and load previously computed Hamiltonians and so
system.target_state = 1
t_g = 20.0
system.target_g = t_g .* [0, 1 ,1, 1, 0, 0]
system.glist = [zeros(151) range(0,50,151) range(0,50,151) zeros(151) zeros(151) zeros(151)]
system4.glist = [1*ones(151) range(0,20,151) 5*ones(151)]
system.numstates = 20
system4.numstates = 20
system3.target_g = [0.0, 10.0, 0.0]
system4.target_g = [1.0, 0.0, 5.0]

println("Computing energy spectrum")
vals = energy_spectrum(system)
vals4 = energy_spectrum(system4)
#gs = system.gcom    # maybe it should be an output of the function, not a parameter of the system
#vals = system.energies # it is stored in a complicated format

xlist = collect(-5.0:0.1:5.0)

println("computing density")
dens = density_profile(system, x = xlist)
dens3 = density_profile(system3, x = xlist)
dens4 = density_profile(system4, x = xlist)

println("computing OBDM eigvals")
obdm_eig = One_Body_Density_Matrix_eigvals(system)
obdm_eig3 = One_Body_Density_Matrix_eigvals(system3)
obdm_eig4 = One_Body_Density_Matrix_eigvals(system4)
println("computing the spatial OBDM")
obdm_s = One_Body_Density_Matrix_spatial(system, x = xlist)
obdm_s3 = One_Body_Density_Matrix_spatial(system3, x = xlist)
obdm_s4 = One_Body_Density_Matrix_spatial(system4, x = xlist)


println("computing the viriral")
virial_e = virial_energies(system)



plot(title = "Energy spectrum",
     xlabel = "g",
     ylabel = "Energy",
     legend = nothing,
     xlims = (0, 50))
for i in 1:system.numstates
    plot!(vals.g[i,:,2], vals.energy[:,i], label = "State $i")
end
display(plot!())

plot(title = "Energy spectrum "*L"g_A=1\; g_B=5",
     xlabel = L"g_{AB}",
     ylabel = "Energy",
     legend = nothing,
     xlims = (0, 20))
for i in 1:system.numstates
    plot!(vals4.g[i,:,2], vals4.energy[:,i], label = "State $i")
end
display(plot!())

#  Now, I store the xlist in the density variable. Mayve we can avoid that. 
#    If we use it as a structure (like dens.x , dens.rho, dens.rhototal ...)
plot(title = "density profile",
     xlabel = "x",
     ylabel = "Density",
     legend = :topright,
     xlims = (-5, 5),
     ylims = (0, 1.1))
for i in 1:3
    plot!(dens.x, dens.dens[i, :], label = "State $i")
end
plot!(dens.x,dens.denstotal,label="Total")
display(plot!())

plot(title = "density profile",
     xlabel = "x",
     ylabel = "Density",
     legend = :topright,
     xlims = (-5, 5),
     ylims = (0, 1.1))
for i in 1:2
    plot!(dens3.x, dens3.dens[i, :], label = "State $i")
end
plot!(dens3.x,dens3.denstotal,label="Total")
display(plot!())

plot(title = "density profile",
     xlabel = "x",
     ylabel = "Density",
     legend = :topright,
     xlims = (-5, 5),
     ylims = (0, 1.1))
for i in 1:2
    plot!(dens4.x, dens4.dens[i, :], label = "State $i")
end
plot!(dens4.x,dens4.denstotal,label="Total")
display(plot!())

# not so much to say here
plot(title = "One Body Density Matrix Eigenvalues",
     xlabel = "Eigenvalue index",
     ylabel = "Eigenvalue",
     legend = nothing,
     ylims = (0, 1))
scatter!(1:length(obdm_eig), obdm_eig, label = "Eigenvalues")
display(plot!())
plot(title = "One Body Density Matrix Eigenvalues",
     xlabel = "Eigenvalue index",
     ylabel = "Eigenvalue",
     legend = nothing,
     ylims = (0, 1))
scatter!(1:length(obdm_eig3), obdm_eig3, label = "Eigenvalues")
display(plot!())
plot(title = "One Body Density Matrix Eigenvalues",
     xlabel = "Eigenvalue index",
     ylabel = "Eigenvalue",
     legend = nothing,
     ylims = (0, 1))
scatter!(1:length(obdm_eig4), obdm_eig4, label = "Eigenvalues")
display(plot!())

#   Here I store the grid in the first index, and the sum on the last one. 
#     Maybe we can avoid that too.
plot(title = "One Body Density Matrix, particle 1",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,)
heatmap!(xlist, xlist, obdm_s[2,:,:], color=:blues, title="OBDM for state 1")
display(plot!())

plot(title = "One Body Density Matrix, particle 2",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,
     )
heatmap!(xlist, xlist, obdm_s[3,:,:], color=:blues, title="OBDM for state 2")
display(plot!())

plot(title = "One Body Density Matrix, particle 1",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,)
heatmap!(xlist, xlist, obdm_s3[2,:,:], color=:blues, title="OBDM for state 1")
display(plot!())

plot(title = "One Body Density Matrix, particle 2",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,
     )
heatmap!(xlist, xlist, obdm_s3[3,:,:], color=:blues, title="OBDM for state 2")
display(plot!())

plot(title = "One Body Density Matrix, particle 1",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,)
heatmap!(xlist, xlist, obdm_s4[2,:,:], color=:blues, title="OBDM for state 1")
display(plot!())

plot(title = "One Body Density Matrix, particle 2",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,
     )
heatmap!(xlist, xlist, obdm_s4[3,:,:], color=:blues, title="OBDM for state 2")
display(plot!())

nothing
