include("../Develop/Hamiltonian_creation.jl")
include("../Develop/Diagonalization.jl")
include("../Develop/Basis_creation.jl")
include("../Develop/Observables.jl")

import .Diagonalization_Develop as Diag
import .HamiltonianCreation_Develop as HC
import .BasisCreation_Develop as BC
import .Observables_Develop as Obs

using Plots


particles = [1, 1]
nho = 10
nstates = 5
type = "bose"
base, dic_base = BC.do_basis(particles, nho, type)
H_ho = HC.one_body_Hamiltonian(particles, nho, type)
H_ID = HC.two_body_Hamiltonian(particles, nho, type)
glist = Vector(0:0.1:50)
vals = Diag.energy_spectrum(H_ho, H_ID, glist, nstates, particles, nho)
println(size(vals))

plot(xlim = (-0.5, 10), xlabel = "g", ylabel = "Energy")
for i in 1:5
    plot!(vals[i, 1, :], vals[i, 2, :], label="State $i")
end
display(plot!())

gt = 20.0
state_t = 2
vecs1 = Diag.eigenstate_calculation(H_ho, H_ID, gt, state_t, particles, nho, corrected=true)
vecs2 = Diag.eigenstate_calculation(H_ho, H_ID, gt, state_t, particles, nho, corrected=false)
println(vecs1 .- vecs2)
dens = Obs.density()
println(dens)

nothing