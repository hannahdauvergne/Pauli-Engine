include("Develop/Struct_main.jl")
using .HO_struct
using Plots

particles = [2,1] 
systemt1 = Few_Particles_Hamonic_Oscillator(particles,30, "bose",parity_restriction=true,parity=true)
systemt2 = Few_Particles_Hamonic_Oscillator(particles,30, "bose",parity_restriction=true,parity=false)
systemt3 = Few_Particles_Hamonic_Oscillator(particles,30, "bose")

l1 = length(systemt1.system.basis)
l2 = length(systemt2.system.basis)
l3 = length(systemt3.system.basis)
println(l1)
println(l2)
println(l3)
println(l3-l2-l1)


systemt1.glist = [ range(0,20,101)  range(0,20,101)  range(0,20,101)]
systemt2.glist = [ range(0,20,101)  range(0,20,101)  range(0,20,101)]
systemt3.glist = [ range(0,20,101)  range(0,20,101)  range(0,20,101)]
systemt1.numstates = 20
systemt2.numstates = 20
systemt3.numstates = 40

energy1 = energy_spectrum(systemt1)
energy2 = energy_spectrum(systemt2)
energy3 = energy_spectrum(systemt3)


plot(xlims = (0,20),ylims = (0,10))

for i in 1:20
    plot!(energy1.g[i,:,2], energy1.energy[:,i], label = "",c=:blue,lw=3)
end
for i in 1:20
    plot!(energy2.g[i,:,2], energy2.energy[:,i], label = "",c=:green,lw=3)
end
for i in 1:40
    plot!(energy3.g[i,:,2], energy3.energy[:,i], label = "",c=:red,ls=:dash,lw=2)
end
display(plot!())

