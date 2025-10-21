############################################   #  N_orbitals = 60
########## function                   time       allocations
##########  mint                      340 ms      3         99 MiB
##########  mint_sparse               5.9 s       28        10 MiB
##########  mint_sparse_v2            356 ms      56        24 MiB
##########  mint_recursive_dict       378 ms      2384558   351 MiB
##########  mint_recursive_self_if    304 ms      2344777   313 MiB
##########  mint_recursive_self       288 ms      2384515   317 MiB
##########  mint_recursive_self_v2    233 ms      1192264   313 MiB  # for me it's the way to go: faster and more precise (seems more allocations than the previous method)
##########  mint_recursive_self_v3    1.5 s       244725682 2.1 GiB
##########  mint_sparse_recursive_dict   6.2 s    2384583   258 MiB
##########  mint_sparse_recursive_dict_v2  362 ms 2384611   270 MiB
##########  mint_sparse_recursive_self_if  6.2 s  2344803   225 MiB
##########  mint_sparse_recursive_self   6.2 s    2384540   228 MiB
###########################################


using Pkg
Pkg.activate("./")


using BenchmarkTools
using SparseArrays
using IterTools
using Combinatorics: permutations

include("../Develop/MintCreation.jl")
using .MintCreation_Develop 


N_orbitals = 60

println("")
println("Compute mint")
@btime compute_mint($N_orbitals)
println("Compute mint sparse")
@btime compute_mint_sparse($N_orbitals)
println("Compute mint sparse v2")
@btime compute_mint_sparse_v2($N_orbitals)
println("Compute mint recursive dict")
@btime compute_mint_recursive_dict($N_orbitals)
println("Compute mint recursive self if")
@btime compute_mint_recursive_self_if($N_orbitals)
println("Compute mint recursive self")
@btime compute_mint_recursive_self($N_orbitals)
println("Compute mint recursive self v2")
@btime compute_mint_recursive_self_v2($N_orbitals)
println("Compute mint recursive self v3")
@btime compute_mint_recursive_self_v3($N_orbitals)
println("Compute mint sparse recursive dict")
@btime compute_mint_sparse_recursive_dict($N_orbitals)
println("Compute mint sparse recursive dict v2")
@btime compute_mint_sparse_recursive_dict_v2($N_orbitals)
println("Compute mint sparse recursive self if")
@btime compute_mint_sparse_recursive_self_if($N_orbitals)
println("Compute mint sparse recursive self")
@btime compute_mint_sparse_recursive_self($N_orbitals)


mint = compute_mint(N_orbitals)
mint_sparse = compute_mint_sparse_v2(N_orbitals)
mint_rec = compute_mint_recursive_dict(N_orbitals)
mint_rec_self_if = compute_mint_recursive_self_if(N_orbitals)
mint_rec_self = compute_mint_recursive_self(N_orbitals)
mint_rec_self = compute_mint_recursive_self_v2(N_orbitals)
mint_sparse_rec = compute_mint_sparse_recursive_dict_v2(N_orbitals)
mint_sparse_rec_self_if = compute_mint_sparse_recursive_self_if(N_orbitals)
mint_sparse_rec_self = compute_mint_sparse_recursive_self(N_orbitals)

println("")
println("Max error mint_rec: ")
println(maximum(abs.(mint .- mint_rec)))
println(argmax(abs.(mint .- mint_rec)))
println(mint[argmax(abs.(mint .- mint_rec))])
println(mint_rec[argmax(abs.(mint .- mint_rec))])
println()
println("Max error mint_rec_self_if: ")
println(maximum(abs.(mint .- mint_rec_self_if)))
println(argmax(abs.(mint .- mint_rec_self_if)))
println(mint[argmax(abs.(mint .- mint_rec_self_if))])
println(mint_rec_self_if[argmax(abs.(mint .- mint_rec_self_if))])
println()
println("Max error mint_rec_self: ")
println(maximum(abs.(mint .- mint_rec_self)))
println(argmax(abs.(mint .- mint_rec_self)))
println(mint[argmax(abs.(mint .- mint_rec_self))])
println(mint_rec_self[argmax(abs.(mint .- mint_rec_self))])
println()
println("Max error mint_sparse_rec: ")
println(maximum(abs.(mint_sparse .- mint_sparse_rec)))
println(argmax(abs.(mint_sparse .- mint_sparse_rec)))
println(mint_sparse[argmax(abs.(mint_sparse .- mint_sparse_rec))])
println(mint_sparse_rec[argmax(abs.(mint_sparse .- mint_sparse_rec))])
println("Max error mint_sparse_rec_self_if: ")
println(maximum(abs.(mint_sparse .- mint_sparse_rec_self_if)))
println(argmax(abs.(mint_sparse .- mint_sparse_rec_self_if)))
println(mint_sparse[argmax(abs.(mint_sparse .- mint_sparse_rec_self_if))])
println(mint_sparse_rec_self_if[argmax(abs.(mint_sparse .- mint_sparse_rec_self_if))])
println("Max error mint_sparse_rec_self: ")
println(maximum(abs.(mint_sparse .- mint_sparse_rec_self)))
println(argmax(abs.(mint_sparse .- mint_sparse_rec_self)))
println(mint_sparse[argmax(abs.(mint_sparse .- mint_sparse_rec_self))])
println(mint_sparse_rec_self[argmax(abs.(mint_sparse .- mint_sparse_rec_self))])

