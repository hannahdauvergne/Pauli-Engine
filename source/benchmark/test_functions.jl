using Pkg
Pkg.activate("./")

include("../Backup/Basis_creation.jl")
include("../Backup/Hamiltonian_creation.jl")
include("../Develop/Basis_creation.jl")
include("../Develop/Hamiltonian_creation.jl")
include("../Develop/MintCreation.jl")
include("../Develop/Diagonalization.jl")

import .BasisCreation as Bc
import .BasisCreation_Develop as BcD
import .HamiltonianCreation as Hc
import .HamiltonianCreation_Develop as HcD
import .MintCreation_Develop as MiCD
import .Diagonalization_Develop as Diag

using Plots
using SpecialFunctions
using BenchmarkTools
using ArnoldiMethod

t0=time()
particles = [1,1]
nho = 10
nstates = 2
type = "bose"
base , dic_base  = BcD.do_basis(particles, nho, type)
H_ho = HcD.one_body_Hamiltonian([1], 2, type)
H_ID = HcD.two_body_Hamiltonian([1], 2, type)
Diag.eigenvalues(H_ho ,H_ID, [1.0], 1)
t2=time()
#H_I = Hc.two_body_Hamiltonian(particles,nho,ptype=type)
H_ho = HcD.one_body_Hamiltonian(particles, nho, type)
H_ID = HcD.two_body_Hamiltonian(particles, nho, type)
#@time HcD.one_body_Hamiltonian(particles, nho, type)
#@time HcD.two_body_Hamiltonian(particles, nho, type)
t1=time()
glist = Vector(-5:0.1:10)
vals = Diag.eigenvalues(H_ho,H_ID, glist,nstates)
vals_p = Diag.eigenvalues_params(H_ho,H_ID, glist,nstates)
vals_p_i = Diag.eigenvalues_params_improved(H_ho,H_ID, glist,nstates)
vals_i = Diag.eigenvalues_improved(H_ho,H_ID, glist,nstates)
vals_a = Diag.eig_arnoldi(H_ho,H_ID, glist,nstates)
vals_a_i = Diag.eig_arnoldi_improved(H_ho,H_ID, glist,nstates)
vals_a_i2 = Diag.eig_arnoldi_improved_v2(H_ho,H_ID, glist,nstates)
c_vals = Diag.correction_eigenvalues(vals, nho, particles, type="fast")
c_vals_a = Diag.correction_eigenvalues(vals_a, nho, particles, type="fast")
c_vals_a_i = Diag.correction_eigenvalues(vals_a_i, nho, particles, type="fast")
c_vals_a_i2 = Diag.correction_eigenvalues(vals_a_i2, nho, particles, type="fast")
println("Time Arpack")
@btime Diag.eigenvalues(H_ho,H_ID, glist,nstates)
println("Time Arpack Improved")
@btime Diag.eigenvalues_improved(H_ho,H_ID, glist,nstates)
println("Time Arpack with parameters")
@btime Diag.eigenvalues_params(H_ho,H_ID, glist,nstates)
println("Time Arpack with parameters improved")
@btime Diag.eigenvalues_params_improved(H_ho,H_ID, glist,nstates)
println("Time Arnoldi")
@btime Diag.eig_arnoldi(H_ho,H_ID, glist,nstates)
println("Time Arnoldi Improved")
@btime Diag.eig_arnoldi_improved(H_ho,H_ID, glist,nstates)
println("Time Arnoldi Improved v2")
@btime Diag.eig_arnoldi_improved_v2(H_ho,H_ID, glist,nstates)


Eex = -2:0.05:1.95
nu = (Eex .-1)./2
gex = -sqrt(8) .* gamma.(1/2 .-nu) ./ gamma.(-nu)

plot(size=(1600, 1200), legend=false)
for i in 1:nstates
    plot!(glist, vals[:,i+1],c=:blue)
    plot!(c_vals[i,1,:] , c_vals[i,2,:], c=:red)
    plot!(c_vals_a[i,1,:] , c_vals_a[i,2,:], c=:green)
    plot!(c_vals_a_i[i,1,:] , c_vals_a_i[i,2,:], c=:orange)
    plot!(c_vals_a_i2[i,1,:] , c_vals_a_i2[i,2,:], c=:lightblue) 
end
plot!(xlims=(-5,10))
scatter!(gex, Eex, markersize=5, c=:black)
plot!(xlabel="g",ylabel="E")
display(plot!())
println(time()-t0)
println(time()-t1)
println(time()-t2)
println(t2-t0)

println(maximum(abs.(c_vals-c_vals_a)))
println(maximum(abs.(c_vals-c_vals_a_i)))