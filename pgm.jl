# Proximal Gradient Method

#=
PGM algorithm.
The module includes the following components:

PGM:       type that holds the PGM variables

setupPGM:  function that initializes the algorithmic variables

AlgoPGM:   function that runs the algorithm


@author : giorgos stathopoulos

@date : 2017-06-20
=#

module ProximalGradient

include("agentCB.jl")
include("agentBESS.jl")
include("master.jl")
using AgentsCB, AgentsBESS, Masters, JuMP, Gurobi

export PGM

type PGM
    u::Array{Any}{1}
    y::Array{Any}{1}
    x::Array{Any}{1}
    v::Array{Any}{1}
    P_CB::Array{Float64, 1}
    P_CB_prev::Array{Float64, 1}
    P_BESS::Array{Float64, 1}
    P_BESS_prev::Array{Float64, 1}
    RES::Array{Float64, 1}
end

function setupPGM(AgentCB::Array{Any,1}, N::Int64, T::Int64)

    # initialize
    pgm_u = Array{Any}(N)
    pgm_y = Array{Any}(N)
    pgm_x = Array{Any}(N)
    pgm_v = Array{Any}(N)

    for i = 1:N
        pgm_u[i] = zeros(AgentCB[i].system.Nu,T)
        pgm_y[i] = zeros(AgentCB[i].system.Ny,T)
        pgm_x[i] = zeros(AgentCB[i].system.Nx,T)
        pgm_v[i] = zeros(T)  # agent's consumption
    end
    pgm_P_CB = zeros(N*T)  # agents' consumptions
    pgm_P_CB_prev = ones(N*T)
    pgm_P_BESS = zeros(T)  # battery's consumption
    pgm_P_BESS_prev = ones(T)
    pgm_RES = [0.0]

    return pgm_u, pgm_y, pgm_x, pgm_v, pgm_P_CB, pgm_P_CB_prev, pgm_P_BESS, pgm_P_BESS_prev, pgm_RES
end

function AlgoPGM(pgm::PGM, N::Int64, T::Int64, Tend::Float64,
    SOC0::Float64, ExecTimes::Dict{Int64,Float64},
    optimizationModelCB::Array{Any,1}, optimizationModelBESS::ProxBESS,
    agentBESS::BESS, agentCB::Array{Any,1}, master::Master, p_BESS_opt::Array{Float64,1}, v_p_CB_opt::Array{Float64,1},
    η::Float64, β::Float64, δ::Float64, α::Float64)

    tempCB   = zeros(N*T)
    tempBESS = zeros(T)

    TimeNextUpdate = collect(values(ExecTimes))[indmax(collect(values(ExecTimes)))]

    k = 1
    T1 = [0.0]
    count = 0

    # main loop
    while T1[end] < Tend

        # Master computes residual
        res = computeResidual(master, α, pgm.P_BESS, pgm.P_CB)

        # BESS proximal step
        pgm.P_BESS_prev = pgm.P_BESS
        point = computeProxPointBESS(agentBESS, δ, β, pgm.P_BESS, pgm.P_BESS_prev, res)
        pgm.P_BESS = getvalue(solveProxBESS(optimizationModelBESS, point, sparse(agentBESS.data.Γ)))

        # CB proximal step
        for i = 1:N
            point = computeProxPointCB(agentCB[i], δ, β, pgm.P_CB[(i-1)*T+1:i*T], pgm.P_CB_prev[(i-1)*T+1:i*T], res)
            pgm.u[i], pgm.y[i], pgm.x[i], pgm.v[i] = solveProxCB(optimizationModelCB[i], point, sparse(agentCB[i].data.Γ))
            tempCB[(i-1)*T+1:i*T] = getvalue(pgm.v[i])
        end

        # update consumption vector
        pgm.P_CB_prev = pgm.P_CB
        pgm.P_CB = tempCB

        # residual
        push!(pgm.RES, norm(pgm.P_BESS - p_BESS_opt)/norm(p_BESS_opt) + norm(pgm.P_CB-v_p_CB_opt)/norm(v_p_CB_opt))

        # save time values
        push!(T1, T1[end]+TimeNextUpdate)
        if mod(k, 10) == 0
            display("Time is: $(T1[end])")
            display("Residual: $(pgm.RES[end])")
        end
        k += 1

    end

    return pgm, T1

end

end
