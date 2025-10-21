include("Develop/Struct_main.jl")
using .HO_struct
using Plots

particles = [2] 
system = Few_Particles_Hamonic_Oscillator(particles,15, "bose")

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

ener1,dens1,norm1  = time_evolution(system,gom1,ts,type="breath",omegas=oms1)
ener2,dens2,norm2  = time_evolution(system,gom2,ts,type="breath",omegas=oms2)
ener3,dens3,norm3  = time_evolution(system,gom3,ts,type="breath",omegas=oms3)
ener4,dens4,norm4  = time_evolution(system,gom4,ts,type="breath",omegas=oms4)


tsf = range(0.0,400000.0,4*np)
ener = [ener1; ener2; ener3; ener4]
dens = [dens1; dens2; dens3; dens4]
norm = [norm1; norm2; norm3; norm4]


plot(tsf,ener)
display(plot!())

plot(tsf,norm)
display(plot!())

heatmap(dens)
display(plot!())


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


tsf = range(0.0,400000.0,4*np)
ener = [ener1; ener2; ener3; ener4]
dens = [dens1; dens2; dens3; dens4]
norm = [norm1; norm2; norm3; norm4]


plot(tsf,ener)
display(plot!())

plot(tsf,norm)
display(plot!())

heatmap(dens)
display(plot!())

#xs = collect(range(-5,5,1001))
#plot(xs,dens[end,:])
#plot!(xs,exp.(-xs .^2)*2/sqrt(pi))
#display(plot!())



