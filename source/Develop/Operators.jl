"""
Implementation of the creation and anihilation operator for both fermions and bosons.

It assumes that the states take the form of a vector ... (explain)
"""
module Operators_Develop
export fop_an,fop_cr,bop_an,bop_cr

function fop_an(state::Vector{Int},indx::Int)
    if indx in state
        phase = (-1)^(findfirst(==(indx), state)-1)
        state[state .== indx] .= 0
        state = sort(state,rev=true)
        return state,phase
    else
        return state,0
    end
end

function fop_cr(state::Vector{Int},indx::Int)
    if indx in state
        return state,0
    else
        state[end] = indx
        state = sort(state,rev=true)
        phase = (-1)^(findfirst(==(indx), state)+1)
        return state,phase
    end
end


function bop_an(state::Vector{Int},indx::Int)
    if indx in state
        count_indx = count(==(indx), state)
        state[findfirst(==(indx), state)] = 0
        state = sort(state,rev=true)
        return state,sqrt(count_indx)
    else
        return state,0
    end
end

function bop_cr(state::Vector{Int},indx::Int)
    state[end] = indx
    state = sort(state,rev=true)
    count_indx = count(==(indx), state)
    return state,sqrt(count_indx)
end

end

