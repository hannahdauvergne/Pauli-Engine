import Pkg 
Pkg.activate("./")

using Test
using Plots

# Import the new versions
include("../Develop/Hamiltonian_creation.jl")
import .HamiltonianCreation_Develop as HC
include("../Develop/Diagonalization.jl")
import .Diagonalization_Develop as Diag

using BenchmarkTools

particles = [1, 1, 1]
nho = 10
nstates = 10
type = "bose"



gl = collect(0:0.1:0.4)
glist = zeros(5,6)
glist[:,2] = gl
glist[:,3] = gl
glist[:,6] = gl

println("Hamiltonian creation")
H_ho = HC.one_body_Hamiltonian(particles, nho, type)
H_ID = HC.two_body_Hamiltonian(particles, nho, type)

#   g in a matrix format only works for the sorted arnoldi case (for now)
gl = collect(0:0.1:15)
glist = zeros(size(gl)[1],6)
glist[:,2] = gl
glist[:,3] = gl
glist[:,6] = gl
valsreq = Diag.energy_spectrum(H_ho,H_ID,glist,nstates,particles,nho,request=true)
valsreq2 = Diag.energy_spectrum(H_ho,H_ID,glist,nstates,particles,nho,request=true,fast=false)

gl = collect(0:0.1:10)
glist = zeros(size(gl)[1],6)
glist[:,2] = gl
glist[:,3] = gl
glist[:,6] = gl
valsnreq = Diag.energy_spectrum(H_ho,H_ID,glist,nstates,particles,nho,request=false)

plot(xlabel="g",ylabel="E",xrange=(0,15),size=(1000,1000))
for i in 1:10
    plot!(valsreq[:,2],valsreq[:,6+i],color=:red,label=:none,lw=3)
end
for i in 1:10
    plot!(valsreq2[:,2],valsreq2[:,6+i],color=:blue,label=:none,lw=3)
end
for i in 1:10
    plot!(valsnreq[i,2,:],valsnreq[i,7,:],c=:lightgreen,ls=:dash,label=:none,lw=3)
end
display(plot!())

println(valsreq[end,7:end])
println(valsreq2[end,7:end])


@btime Diag.energy_spectrum(H_ho,H_ID,glist,nstates,particles,nho,request=true)            # controls the range but the results are inacurate (ground state is almost perfect corrected)
@btime Diag.energy_spectrum(H_ho,H_ID,glist,nstates,particles,nho,request=true,fast=false) # The result is nice, but it is so slow (x40 than others)
@btime Diag.energy_spectrum(H_ho,H_ID,glist,nstates,particles,nho,request=false)           # The output range does not match with the input one
nothing