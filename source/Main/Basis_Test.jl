# Basis_Test.jl
include("Basis.jl")

println("======================================")
println("   Testing basis_creation (Dict basis)")
println("======================================")

"""
    run_test(particles, lvls, particletype)

Run a basis_creation test and print results.
"""
function run_test(particles, lvls, particletype)
    println("\n--------------------------------------")
    println("Test: particles = $particles, lvls = $lvls, type = $particletype")
    println("--------------------------------------")

    basis = basis_creation(particles, lvls, particletype)

    println("\nBasis dictionary (state → index):")
    for (k, v) in basis
        println("$k  →  $v")
    end

    println("\nNumber of basis states: ", length(basis))

    println("\nIndex → state check (sorted by index):")
    inv = Dict(v => k for (k, v) in basis)
    for i in 1:length(inv)
        println("Index $i → ", inv[i])
    end
end


###############################
#         TEST CASES
###############################

# 1) One boson, 3 levels
run_test([1], 3, "bose")

# 2) Two bosons, 3 levels
run_test([2], 3, "bose")

# 3) One fermion, 3 levels
run_test([1], 3, "fermi")

# 4) Two fermions, 4 levels
run_test([2], 4, "fermi")

# 5) Three bosons, 3 levels
run_test([3], 3, "bose")

# 6) Three fermions, 3 levels
run_test([3], 3, "fermi")


println("\nDone.")
