include("Main/Main.jl")
using .Harmonic_Oscillator_solver
using Plots

particles = [2] 
system = Few_Particles_Hamonic_Oscillator(particles,10, "bose")

np = 101
gs = range(0.0,15.0,np)
ts = range(0.0,10000.0,np)

oms1 = range(0.0,0.0,np)  # omega_new^2 - 1
gom1 = range(0.0,15.0,np) 

oms2 = range(0.0,8.0,np)  # omega_new^2 - 1
gom2 = range(15.0,15.0,np) 

oms3 = range(8.0,8.0,np)  # omega_new^2 - 1
gom3 = range(15.0,0.0,np) 

oms4 = range(8.0,0.0,np)  # omega_new^2 - 1
gom4 = range(0.0,0.0,np) 

ts4 = range(0.0,40000.0,401)
gom = [gom1; gom2; gom3; gom4]
oms = [oms1; oms2; oms3; oms4]

ener1,dens1,norm1 = time_evolution(system,gom1,ts,type="breath",omegas=oms1)
ener2,dens2,norm2 = time_evolution(system,gom2,ts,type="breath",omegas=oms2)
ener3,dens3,norm3 = time_evolution(system,gom3,ts,type="breath",omegas=oms3)
ener4,dens4,norm4 = time_evolution(system,gom4,ts,type="breath",omegas=oms4)

enert,denst,normt = time_evolution(system,gom,ts4,type="breath",omegas=oms)

tsf = range(0.0,40000.0,4*np)
ener = [ener1; ener2; ener3; ener4]
dens = [dens1; dens2; dens3; dens4]
norm = [norm1; norm2; norm3; norm4]


plot(tsf,ener,label="")
plot!(ts4,enert,label="")
plot!(xlabel="t",ylabel="E")
display(plot!())
savefig("Figs/energy_vs_t.pdf")

plot(tsf,norm,label="")
plot!(ts4,normt,label="")
plot!(xlabel="t",ylabel="Norm")
display(plot!())
savefig("Figs/norm_vs_t.pdf")

xs = collect(range(-5,5,1001))
heatmap(xs,tsf,dens,clim=(0,2.0))
plot!(xlabel="x",ylabel="t")
display(plot!())
savefig("Figs/density_vs_t.pdf")

heatmap(xs,ts4,denst,clim=(0.0,2.0))
plot!(xlabel="x",ylabel="t")
display(plot!())
savefig("Figs/density_vs_t_continuous.pdf")


oms1 = range(0.0,8.0,np)  # omega_new^2 - 1
gom1 = range(0.0,0.0,np) 

oms2 = range(8.0,8.0,np)  # omega_new^2 - 1
gom2 = range(0.0,15.0,np) 

oms3 = range(8.0,0.0,np)  # omega_new^2 - 1
gom3 = range(15.0,15.0,np) 

oms4 = range(0.0,0.0,np)  # omega_new^2 - 1
gom4 = range(15.0,0.0,np) 

ener1,dens1,norm1  = time_evolution(system,gom1,ts,type="breath",omegas=oms1)
ener2,dens2,norm2  = time_evolution(system,gom2,ts,type="breath",omegas=oms2)
ener3,dens3,norm3  = time_evolution(system,gom3,ts,type="breath",omegas=oms3)
ener4,dens4,norm4  = time_evolution(system,gom4,ts,type="breath",omegas=oms4)


tsfv2 = range(0.0,40000.0,4*np)
enerv2 = [ener1; ener2; ener3; ener4]
densv2 = [dens1; dens2; dens3; dens4]
normv2 = [norm1; norm2; norm3; norm4]


plot(tsfv2,enerv2,label="")
plot!(xlabel="t",ylabel="E")
display(plot!())
savefig("Figs/energy_vs_t_reverse.pdf")

plot(tsfv2,normv2,label="")
plot!(xlabel="t",ylabel="Norm")
display(plot!())
savefig("Figs/norm_vs_t_reverse.pdf")

heatmap(xs,tsfv2,densv2)
plot!(xlabel="x",ylabel="t")
display(plot!())
savefig("Figs/density_vs_t_reverse.pdf")

#xs = collect(range(-5,5,1001))
#plot(xs,dens[end,:])
#plot!(xs,exp.(-xs .^2)*2/sqrt(pi))
#display(plot!())



