import Pkg 
Pkg.activate("./")

# Import the new versions
include("../Develop/Hamiltonian_creation.jl")
import .HamiltonianCreation_Develop as HCD
include("../Backup/Hamiltonian_creation.jl")
import .HamiltonianCreation as HC

using BenchmarkTools

particles = [2]
ptype = "bose"
nho = 40

H1 = HC.one_body_Hamiltonian(particles,nho,ptype)
println("time state loop (1 body)")
@btime HC.one_body_Hamiltonian(particles,nho,ptype)

H1T = HCD.one_body_Hamiltonian(particles,nho,ptype)
println("time dictionary loop (1 body)")
@btime HCD.one_body_Hamiltonian(particles,nho,ptype)

println(H1 == H1T)


H2 = HC.two_body_Hamiltonian(particles,nho,ptype)
println("time state loop (2 body)")
@btime HC.two_body_Hamiltonian(particles,nho,ptype)

H2T = HCD.two_body_Hamiltonian(particles,nho,ptype)
println("time dictionary loop (2 body)")
@btime HCD.two_body_Hamiltonian(particles,nho,ptype)

println(H2 == H2T)

