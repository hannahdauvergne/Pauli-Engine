"""
        Load all the modules needed (maybe can be loaded in each other files), as well the different files.
        Export all the functions that the user could need.  
"""

module Harmonic_Oscillator_solver
using SpecialPolynomials, Polynomials
using SpecialFunctions
using SparseArrays
using ArnoldiMethod
using JLD2

include("Model.jl")
include("Basis.jl")
include("Hamiltonian.jl")
include("Diagonalization.jl")
include("Observables.jl")
include("Time_evolution.jl")


export Few_Particles_Hamonic_Oscillator,one_body_Hamiltonian!,two_body_Hamiltonian!,eigenstate!,
        energy_spectrum,density_profile,One_Body_Density_Matrix,One_Body_Density_Matrix_eigvals,
        One_Body_Density_Matrix_spatial,virial_energies,momentum_distribution,pair_correlation,
        pair_correlation_spatial,quench,time_evolution

end
