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

include("model.jl")
include("basis.jl")
include("Hamiltonian.jl")
include("diagonalization.jl")
include("observables.jl")



export Few_Particles_Hamonic_Oscillator,one_body_Hamiltonian!,two_body_Hamiltonian!,eigenstate!,
        energy_spectrum,density_profile,One_Body_Density_Matrix,One_Body_Density_Matrix_eigvals,
        One_Body_Density_Matrix_spatial,virial_energies

end
