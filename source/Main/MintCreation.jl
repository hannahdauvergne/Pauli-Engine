"""
Different functions to compute the matrix elements of the contact interaction.
We use the "compute_mint_recursive_self_v2 as it reported the best perfomance.
"""
module MintCreation

export compute_mint, compute_mint_sparse, compute_mint_sparse_v2, compute_mint_recursive_dict,
    compute_mint_recursive_self_if, compute_mint_recursive_self, compute_mint_recursive_self_v2,
    compute_mint_recursive_self_v3, compute_mint_sparse_recursive_dict,
    compute_mint_sparse_recursive_dict_v2, compute_mint_sparse_recursive_self_if,
    compute_mint_sparse_recursive_self

using SpecialFunctions
using SparseArrays
using IterTools, Combinatorics


function compute_mint(sites::Int)

    mint = zeros(Float64, sites, sites, sites, sites)
    for i in 1:sites 
        for j in 1:i 
            for k in 1:j 
                for l in ((i+j+k+1)%2)+1:2:k 
                    mint[i,j,k,l] = I_2(i-1,j-1,k-1,l-1)
                end
            end
        end
    end

    return mint

end

function compute_mint_sparse(sites::Int)

    mint = spzeros(sites*sites, sites*sites)
    for i in 1:sites 
        for j in 1:i 
            for k in 1:j 
                for l in ((i+j+k+1)%2)+1:2:k 
                    mint[(i-1)*sites + j, (k-1)*sites + l] = I_2(i-1,j-1,k-1,l-1)
                end
            end
        end
    end

    return mint

end

function compute_mint_sparse_v2(sites::Int)
    rows = Int[]
    columns = Int[]
    values = Float64[]
    for i in 1:sites 
        for j in 1:i 
            for k in 1:j 
                for l in ((i+j+k+1)%2)+1:2:k 
                    push!(rows, (i-1)*sites + j)
                    push!(columns, (k-1)*sites + l)
                    push!(values, I_2(i-1,j-1,k-1,l-1))
                end
            end
        end
    end
    mint = sparse(rows, columns, values, sites*sites, sites*sites)
    return mint

end

function compute_mint_recursive_dict(sites::Int)

    mint = zeros(Float64, sites, sites, sites, sites)
    mint[1,1,1,1] = 1/sqrt(2 * pi)
    W = Dict((1,1,1,1) => 1/sqrt(2 * pi))
    for i in 2:sites 
        for j in 1:i 
            for k in 1:j 
                for l in ((i+j+k+1)%2)+1:2:k 
                    mint[i,j,k,l] = 0.5*(-sqrt((i-2)/(i-1)) * calc_W([i-2,j,k,l], W) +
                    sqrt((j-1)/(i-1)) * calc_W([i-1,j-1,k,l], W) +
                    sqrt((k-1)/(i-1)) * calc_W([i-1,j,k-1,l], W) +
                    sqrt((l-1)/(i-1)) * calc_W([i-1,j,k,l-1], W))
                    index = (i,j,k,l)
                    W[index] = mint[i,j,k,l]
                end
            end
        end
    end

    return mint

end

function compute_mint_recursive_self_if(sites::Int)

    mint = zeros(Float64, sites, sites, sites, sites)
    mint[1,1,1,1] = 1/sqrt(2 * pi)
    for i in 2:sites 
        for j in 1:i 
            for k in 1:j 
                for l in ((i+j+k+1)%2)+1:2:k 
                    if i>2
                        a,b,c,d = sort([i-2,j,k,l],rev=true)
                        mint[i,j,k,l] -= sqrt((i-2)/(i-1)) * mint[a,b,c,d]
                    end
                    if j>1
                        a,b,c,d = sort([i-1,j-1,k,l],rev=true)
                        mint[i,j,k,l] += sqrt((j-1)/(i-1)) * mint[a,b,c,d] 
                    end
                    if k>1
                        a,b,c,d = sort([i-1,j,k-1,l],rev=true)
                        mint[i,j,k,l] += sqrt((k-1)/(i-1)) * mint[a,b,c,d]
                    end
                    if l>1
                        a,b,c,d = sort([i-1,j,k,l-1],rev=true)
                        mint[i,j,k,l] += sqrt((l-1)/(i-1)) * mint[a,b,c,d]
                    end
                    mint[i,j,k,l] *= 0.5
                end
            end
        end
    end

    return mint

end

function compute_mint_recursive_self(sites::Int)

    mint = zeros(Float64, sites, sites, sites, sites)
    mint[1,1,1,1] = 1/sqrt(2 * pi)
    for i in 2:sites 
        for j in 1:i 
            for k in 1:j 
                for l in ((i+j+k+1)%2)+1:2:k 
                    mint[i,j,k,l] = 0.5*(-sqrt((i-2)/(i-1)) * calcmint([i-2,j,k,l], mint) +
                    sqrt((j-1)/(i-1)) * calcmint([i-1,j-1,k,l], mint) +
                    sqrt((k-1)/(i-1)) * calcmint([i-1,j,k-1,l], mint) +
                    sqrt((l-1)/(i-1)) * calcmint([i-1,j,k,l-1], mint))
                end
            end
        end
    end

    return mint

end

function compute_mint_recursive_self_v2(sites::Int)
    mint = zeros(Float64, sites+1, sites+1, sites+1, sites+1)
    mint[2,2,2,2] = 1/sqrt(2 * pi)
    temp_indices = Vector{Int}(undef, 4) # Preallocate for sorting
    sqrt_cache = sqrt.(0:(sites-1)) # Precompute square roots for efficiency
    for i in 3:sites+1 
        sqrt_i_2 = sqrt_cache[i-1]
        sqrt_i_3 = sqrt_cache[i-2]
        for j in 2:i 
            sqrt_j_2 = sqrt_cache[j-1]
            for k in 2:j 
                sqrt_k_2 = sqrt_cache[k-1]
                for l in ((i+j+k)%2)+2:2:k
                    sqrt_l_2 = sqrt_cache[l-1]
                    temp_indices[1:4] = [i-2, j, k, l]
                    a1, b1, c1, d1 = sort!(temp_indices, rev=true)
                    temp_indices[1:4] = [i-1, j-1, k, l]
                    a2, b2, c2, d2 = sort!(temp_indices, rev=true)
                    temp_indices[1:4] = [i-1, j, k-1, l]
                    a3, b3, c3, d3 = sort!(temp_indices, rev=true)
                    temp_indices[1:4] = [i-1, j, k, l-1]
                    a4, b4, c4, d4 = sort!(temp_indices, rev=true)
                    mint[i, j, k, l] = 0.5 * (
                        -sqrt_i_3 / sqrt_i_2 * mint[a1, b1, c1, d1] +
                        sqrt_j_2 / sqrt_i_2 * mint[a2, b2, c2, d2] +
                        sqrt_k_2 / sqrt_i_2 * mint[a3, b3, c3, d3] +
                        sqrt_l_2 / sqrt_i_2 * mint[a4, b4, c4, d4]
                    )
                end
            end
        end
    end
    return mint[2:end,2:end,2:end,2:end]
end

function compute_mint_recursive_self_v2_trunc(sites::Int, pmax :: Int)
    mint = zeros(Float64, sites+1, sites+1, sites+1, sites+1)
    mint[2,2,2,2] = 1/sqrt(2 * pi)
    temp_indices = Vector{Int}(undef, 4) # Preallocate for sorting
    sqrt_cache = sqrt.(0:(sites-1)) # Precompute square roots for efficiency
    for i in 3:sites+1 
        sqrt_i_2 = sqrt_cache[i-1]
        sqrt_i_3 = sqrt_cache[i-2]
        for j in 2:i 
            sqrt_j_2 = sqrt_cache[j-1]
            for k in 2:min(j,2*sites + pmax - i - j + 3) 
                sqrt_k_2 = sqrt_cache[k-1]
                for l in ((i+j+k)%2)+2:2:min(k, 2*sites + pmax - i - j - k + 5)
                    sqrt_l_2 = sqrt_cache[l-1]
                    temp_indices[1:4] = [i-2, j, k, l]
                    a1, b1, c1, d1 = sort!(temp_indices, rev=true)
                    temp_indices[1:4] = [i-1, j-1, k, l]
                    a2, b2, c2, d2 = sort!(temp_indices, rev=true)
                    temp_indices[1:4] = [i-1, j, k-1, l]
                    a3, b3, c3, d3 = sort!(temp_indices, rev=true)
                    temp_indices[1:4] = [i-1, j, k, l-1]
                    a4, b4, c4, d4 = sort!(temp_indices, rev=true)
                    mint[i, j, k, l] = 0.5 * (
                        -sqrt_i_3 / sqrt_i_2 * mint[a1, b1, c1, d1] +
                        sqrt_j_2 / sqrt_i_2 * mint[a2, b2, c2, d2] +
                        sqrt_k_2 / sqrt_i_2 * mint[a3, b3, c3, d3] +
                        sqrt_l_2 / sqrt_i_2 * mint[a4, b4, c4, d4]
                    )
                end
            end
        end
    end
    return mint[2:end,2:end,2:end,2:end]
end

function compute_mint_recursive_self_v3(sites::Int)
    mint = zeros(Float64, sites+1, sites+1, sites+1, sites+1)
    mint[2,2,2,2] = 1/sqrt(2 * pi)
    sqrt_cache = sqrt.(0:(sites-1)) # Precompute square roots for efficiency
    for i in 3:sites+1 
        sqrt_i_2 = sqrt_cache[i-1]
        sqrt_i_3 = sqrt_cache[i-2]
        for j in 2:i 
            sqrt_j_2 = sqrt_cache[j-1]
            for k in 2:j 
                sqrt_k_2 = sqrt_cache[k-1]
                for l in ((i+j+k)%2)+2:2:k
                    sqrt_l_2 = sqrt_cache[l-1]
                    value = 0.5 * (
                        -sqrt_i_3 / sqrt_i_2 * mint[i-2,j,k,l] +
                         sqrt_j_2 / sqrt_i_2 * mint[i-1,j-1,k,l] +
                         sqrt_k_2 / sqrt_i_2 * mint[i-1,j,k-1,l] +
                         sqrt_l_2 / sqrt_i_2 * mint[i-1,j,k,l-1]
                    )
                    # Assign value to all permutations of (i, j, k, l)
                    for (a, b, c, d) in unique(permutations((i, j, k, l)))
                        mint[a, b, c, d] = value
                    end
                end
            end
        end
    end
    return mint[2:end,2:end,2:end,2:end]
end

function calcmint(index::Vector{Int}, mint::Array{Float64,4})
    i,j,k,l = sort(index, rev=true)
    if l > 0
        return mint[i,j,k,l]
    else
        return 0.0
    end
end


# Read the value from the dictionary or return 0 
function calc_W(index::Vector{Int}, W::Dict{NTuple{4,Int}, Float64})
    # All indices are positive
    index = sort(index, rev=true) # Sort indices in descending order
    index = NTuple{4, Int}(index) # Convert the index vector to a tuple of 4 integers
    return get(W, index, 0.0) # Get the value from the dictionary or default to 0.0
end

function compute_mint_sparse_recursive_dict(sites::Int)

    mint = spzeros(sites*sites, sites*sites)
    #mint[1,1,1,1] = 1/sqrt(2 * pi)
    mint[1,1] = 1/sqrt(2 * pi)
    W = Dict((1,1,1,1) => 1/sqrt(2 * pi))
    for i in 2:sites 
        for j in 1:i 
            for k in 1:j 
                for l in ((i+j+k+1)%2)+1:2:k 
                    mint[(i-1)*sites + j,(k-1)*sites + l] = 0.5*(-sqrt((i-2)/(i-1)) * calc_W([i-2,j,k,l], W) +
                    sqrt((j-1)/(i-1)) * calc_W([i-1,j-1,k,l], W) +
                    sqrt((k-1)/(i-1)) * calc_W([i-1,j,k-1,l], W) +
                    sqrt((l-1)/(i-1)) * calc_W([i-1,j,k,l-1], W))
                    index = (i,j,k,l)
                    W[index] = mint[(i-1)*sites + j,(k-1)*sites + l]
                end
            end
        end
    end
    return mint

end

function compute_mint_sparse_recursive_dict_v2(sites::Int)

    rows = Int[1]
    columns = Int[1]
    values = Float64[1/sqrt(2 * pi)]
    W = Dict((1,1,1,1) => 1/sqrt(2 * pi))
    for i in 2:sites 
        for j in 1:i 
            for k in 1:j 
                for l in ((i+j+k+1)%2)+1:2:k 
                    value = 0.5*(-sqrt((i-2)/(i-1)) * calc_W([i-2,j,k,l], W) +
                    sqrt((j-1)/(i-1)) * calc_W([i-1,j-1,k,l], W) +
                    sqrt((k-1)/(i-1)) * calc_W([i-1,j,k-1,l], W) +
                    sqrt((l-1)/(i-1)) * calc_W([i-1,j,k,l-1], W))
                    push!(rows, (i-1)*sites + j)
                    push!(columns, (k-1)*sites + l)
                    push!(values, value)
                    index = (i,j,k,l)
                    W[index] = value
                end
            end
        end
    end
    mint = sparse(rows, columns, values, sites*sites, sites*sites)
    return mint

end

function compute_mint_sparse_recursive_self_if(sites::Int)

    mint = spzeros(sites*sites, sites*sites)
    mint[1,1] = 1/sqrt(2 * pi)
    a,b,c,d = [1,1,1,1]
    for i in 2:sites 
        for j in 1:i 
            for k in 1:j 
                for l in ((i+j+k+1)%2)+1:2:k 
                    if i>2
                        a,b,c,d = sort([i-2,j,k,l],rev=true)
                        mint[(i-1)*sites + j,(k-1)*sites + l] -= sqrt((i-2)/(i-1)) * mint[(a-1)*sites + b,(c-1)*sites + d]
                    end
                    if j>1
                        a,b,c,d = sort([i-1,j-1,k,l],rev=true)
                        mint[(i-1)*sites + j,(k-1)*sites + l] += sqrt((j-1)/(i-1)) * mint[(a-1)*sites + b,(c-1)*sites + d] 
                    end
                    if k>1
                        a,b,c,d = sort([i-1,j,k-1,l],rev=true)
                        mint[(i-1)*sites + j,(k-1)*sites + l] += sqrt((k-1)/(i-1)) * mint[(a-1)*sites + b,(c-1)*sites + d]
                    end
                    if l>1
                        a,b,c,d = sort([i-1,j,k,l-1],rev=true)
                        mint[(i-1)*sites + j,(k-1)*sites + l] += sqrt((l-1)/(i-1)) * mint[(a-1)*sites + b,(c-1)*sites + d]
                    end
                    mint[(i-1)*sites + j,(k-1)*sites + l] *= 0.5
                end
            end
        end
    end

    return mint

end

function compute_mint_sparse_recursive_self(sites::Int)

    mint = spzeros(sites * sites, sites * sites)
    mint[1,1] = 1/sqrt(2 * pi)
    for i in 2:sites 
        for j in 1:i 
            for k in 1:j 
                for l in ((i+j+k+1)%2)+1:2:k 
                    mint[(i-1)*sites + j,(k-1)*sites + l] = 0.5*(-sqrt((i-2)/(i-1)) * calcmint_sparse([i-2,j,k,l], mint, sites) +
                    sqrt((j-1)/(i-1)) * calcmint_sparse([i-1,j-1,k,l], mint, sites) +
                    sqrt((k-1)/(i-1)) * calcmint_sparse([i-1,j,k-1,l], mint, sites) +
                    sqrt((l-1)/(i-1)) * calcmint_sparse([i-1,j,k,l-1], mint, sites))
                end
            end
        end
    end

    return mint

end

function calcmint_sparse(index::Vector{Int}, mint::SparseMatrixCSC{Float64, Int64}, sites::Int)
    i,j,k,l = sort(index, rev=true)
    if l > 0
        return mint[(i-1)*sites+j,(k-1)*sites + l]
    else
        return 0.0
    end
end

function I_2(a,b,c,d)
    integral=0.0
    if ((a+b+c+d)%2 == 0)
      for r in 0:d
        int_r = -0.5*log(2)-2*log(pi)
        int_r += 0.5*(loggamma(c+1)+loggamma(d+1)-loggamma(a+1)-loggamma(b+1))  # gamma(n+1)=n!
        int_r += -loggamma(d-r+1)-loggamma(c-d+r+1)-loggamma(r+1) # gamma(n+1)=n!
        logabs1, sign1 = logabsgamma((a+b-c+d+1)/2-r)
        logabs2, sign2 = logabsgamma((a-b+c-d+1)/2+r)
        logabs3, sign3 = logabsgamma((-a+b+c-d+1)/2+r)
        int_r += logabs1+logabs2+logabs3
        phase_r = sign1*sign2*sign3
        integral += exp(int_r)*phase_r
      end
    end
    return integral
end

end