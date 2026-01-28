##############################
# Test code
##############################
include("Operators.jl")
using .Operators

println("=== FERMIONS ===")

# Fermionic state: modes 5, 3, 1 occupied
f_state = [5, 3, 1, 0, 0]
println("Initial state: ", f_state)

s, p = fop_an(copy(f_state), 3)
println("c₃ |ψ⟩ → state = ", s, ", phase = ", p)

s, p = fop_an(copy(f_state), 2)
println("c₂ |ψ⟩ → state = ", s, ", phase = ", p)

s, p = fop_cr(copy(f_state), 4)
println("c₄† |ψ⟩ → state = ", s, ", phase = ", p)

s, p = fop_cr(copy(f_state), 3)
println("c₃† |ψ⟩ → state = ", s, ", phase = ", p)

println("\n=== BOSONS ===")

# Bosonic state: three bosons in mode 2, one in mode 1
b_state = [4,4,4,3,2,2,0]
println("Initial state: ", b_state)

s, f = bop_an(copy(b_state), 2)
println("a₂ |ψ⟩ → state = ", s, ", factor = ", f)

s, f = bop_cr(copy(b_state), 2)
println("a₂† |ψ⟩ → state = ", s, ", factor = ", f)

println("\n=== FERMION ANTICOMMUTATION CHECK ===")

# |ψ⟩ = c₃† c₁† |0⟩
state = [3, 1, 0, 0]

s1, p1 = fop_an(copy(state), 3)
s1, p1b = fop_an(s1, 1)
A = p1 * p1b

s2, p2 = fop_an(copy(state), 1)
s2, p2b = fop_an(s2, 3)
B = p2 * p2b

println("c₁ c₃ |ψ⟩ phase = ", A)
println("c₃ c₁ |ψ⟩ phase = ", B)
println("Sum: ", A + B)