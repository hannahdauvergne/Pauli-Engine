include("Develop/Struct_main.jl")
using .HO_struct
using Plots
using LaTeXStrings

particles = [2,1] 
systemb = Few_Particles_Hamonic_Oscillator(particles,30, "bose")
systemf = Few_Particles_Hamonic_Oscillator(particles,30, "fermi")

systemb.target_state = 1
systemb.target_g = [ 1000.0 , 1.0 ,0.0]
systemb.glist = [1000.0 .* ones(101) range(0.1,15,101) zeros(101)]
systemb.numstates = 30

systemf.target_state = 1
systemf.target_g = [ 0.0 , 1.0 ,0.0]
systemf.glist = [0.0 .* ones(101) range(0.1,15,101) zeros(101)]
systemf.numstates = 30

xlist = collect(-5.0:0.1:5.0)

println("Computing energy spectrum")
valsb = energy_spectrum(systemb)
valsf = energy_spectrum(systemf)

println("computing density")
densb = density_profile(systemb, x = xlist)
densf = density_profile(systemf, x = xlist)

println("computing OBDM eigvals")
obdm_eigb = One_Body_Density_Matrix_eigvals(systemb)
obdm_eigf = One_Body_Density_Matrix_eigvals(systemf)

println("computing the spatial OBDM")
obdm_sb = One_Body_Density_Matrix_spatial(systemb, x = xlist)
obdm_sf = One_Body_Density_Matrix_spatial(systemf, x = xlist)

println("computing the viriral")
virial_eb = virial_energies(systemb)
virial_ef = virial_energies(systemf)

println("computing the momentum distribution")
momentumb = momentum_distribution(systemb,x = xlist, k = collect(-5.0:0.1:5.0))
momentumf = momentum_distribution(systemf,x = xlist, k = collect(-5.0:0.1:5.0))

println("computing the pair correlation function")
PCFb = pair_correlation(systemb)
PCFf = pair_correlation(systemf)

println("computing the pair correlation spatial")
PCF_sb = pair_correlation_spatial(systemb, x = xlist)
PCF_sf = pair_correlation_spatial(systemf, x = xlist)



plot(title = "Energy spectrum",
     xlabel = "g",
     ylabel = "Energy",
     legend = nothing,
     xlims = (0, 15))
for i in 1:systemb.numstates
    plot!(valsb.g[i,:,2], valsb.energy[:,i], label = "State $i",c=:blue)
    plot!(valsf.g[i,:,2], valsf.energy[:,i], label = "State $i",c=:red,ls=:dash)
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
for i in 1:2
    plot!(densb.x, densb.dens[i, :], label = "",c=:blue)
    plot!(densb.x, densf.dens[i, :], label = "",c=:red,ls=:dash)
end
plot!(densb.x,densb.denstotal,label="",c=:blue)
plot!(densf.x,densf.denstotal,label="",c=:red,ls=:dash)
display(plot!())

# not so much to say here
plot(title = "One Body Density Matrix Eigenvalues",
     xlabel = "Eigenvalue index",
     ylabel = "Eigenvalue",
     legend = nothing,
     ylims = (0, 1))
scatter!(1:length(obdm_eigb), obdm_eigb, label = "Eigenvalues",c=:blue)
scatter!(1:length(obdm_eigf), obdm_eigf, label = "Eigenvalues",c=:red)
display(plot!())

#   Here I store the grid in the first index, and the sum on the last one. 
#     Maybe we can avoid that too.
plot(title = "One Body Density Matrix, particle 1 (b)",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,)
heatmap!(xlist, xlist, obdm_sb[2,:,:], color=:blues, title="OBDM for state 1")
display(plot!())


plot(title = "One Body Density Matrix, particle 1 (f))",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,)
heatmap!(xlist, xlist, obdm_sf[2,:,:], color=:blues, title="OBDM for state 1")
display(plot!())

plot(title = "One Body Density Matrix, particle 2 (b)",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,
     )
heatmap!(xlist, xlist, obdm_sb[3,:,:], color=:blues, title="OBDM for state 2")
display(plot!())

plot(title = "One Body Density Matrix, particle 2 (f)",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,
     )
heatmap!(xlist, xlist, obdm_sf[3,:,:], color=:blues, title="OBDM for state 2")
display(plot!())

plot(title = "Momentum distribution, particle 1",
     xlabel = "k",
     ylabel = "p_k'",
     legend = nothing,
     )
plot!(xlist, momentumb[2,:,:], color=:blue)
plot!(xlist, momentumf[2,:,:], color=:red,ls=:dash)
display(plot!())

plot(title = "Momentum distribution, particle 2",
     xlabel = "k",
     ylabel = "p_k'",
     legend = nothing,
     )
plot!(xlist, momentumb[3,:,:], color=:blue)
plot!(xlist, momentumf[3,:,:], color=:red,ls=:dash)
display(plot!())

plot(title = "Pair correlation 1-1 (b)",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,
     )
heatmap!(xlist, xlist, PCF_sb[1,1,:,:], color=:blues)
display(plot!())

plot(title = "Pair correlation 1-1 (f)",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,
     )
heatmap!(xlist, xlist, PCF_sf[1,1,:,:], color=:blues)
display(plot!())

plot(title = "Pair correlation 1-2 (b)",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,
     )
heatmap!(xlist, xlist, PCF_sb[1,2,:,:], color=:blues)
display(plot!())

plot(title = "Pair correlation 1-2 (f)",
     xlabel = "x",
     ylabel = "x'",
     legend = nothing,
     )
heatmap!(xlist, xlist, PCF_sf[1,2,:,:], color=:blues)
display(plot!())


nothing


