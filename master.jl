# Master data

#=
Sets up the master node.
The module includes the following components:

Master: constructor

setupMaster: type that includes the optimization data

computeResidual: gathers power profiles from agents and
                 computes difference to tracked signal

@author : giorgos stathopoulos

@date : 2017-06-20
=#

module Masters

using DataFrames

export Master, computeResidual

type Master
    N::Int64                   # No of agents
    S::Array{Float64, 2}       # selection matrix
    E::Array{Float64, 2}       # summation matrix
    REF::Array{Float64, 1}     # reference signal (aggregate deviation from baseline)
    ϵ_unc::Array{Float64, 1}   # error estimate that needs to be eliminated
end

function setupMaster(S::Array{Any}, E::Array{Any}, REF::Array{Float64, 1}, N::Int64,  p_CB_tot::Array{Float64, 1})

    # gradient terms; Hessian, selection matrix and summation matrix
    S   = collect(convert(Array,readtable(string("data/SS", ".dat"), separator = ',', header = false)))
    E   = collect(convert(Array,readtable(string("data/EE", ".dat"), separator = ',', header = false)))

    # signal to compensate for
    ϵ_unc = REF + p_CB_tot

    master = Master(N, S, E, REF, ϵ_unc)

end

function computeResidual(master::Master, α::Float64, p_BESS::Array{Float64, 1}, p_CB::Array{Float64, 1})

    # solve problem
    res = α*(p_BESS + master.S*master.E*p_CB - master.ϵ_unc)

    return vec(res)
end


# function computeDual(master::Master, λ::Array{Float64, 1}, p_BESS::Array{Float64, 1}, p_CB::Array{Float64, 1}, ρ::Float64)

#     # compute multiplier
#     z = λ + ρ * (p_BESS - master.S*master.E*p_CB - master.ϵ_unc)

#     return z
# end

end
