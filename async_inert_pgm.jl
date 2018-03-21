# Asynchronous (Inertial) Proximal Gradient Method

#=
AsInPGM algorithm.
The module includes the following components:

AsInPGM:       type that holds the AsInPGM variables

setupAsInPGM:  function that initializes the algorithmic variables

AlgoAsInPGM:   function that runs the algorithm


@author : giorgos stathopoulos

@date : 2017-06-20
=#

module AsynchronousProximalGradient

include("agentCB.jl")
include("agentBESS.jl")
include("master.jl")
using AgentsCB, AgentsBESS, Masters, JuMP, Gurobi

export AsInPGM

type AsInPGM
    u::Array{Any}{1}
    y::Array{Any}{1}
    x::Array{Any}{1}
    v::Array{Any}{1}
    P_CB::Array{Float64, 2}
    P_CB_prev::Array{Float64, 2}
    P_BESS::Array{Float64, 1}
    P_BESS_prev::Array{Float64, 1}
    RES::Array{Float64, 1}
    P_CB_local::Array{Float64, 2}
    P_BESS_local::Array{Float64, 2}
    P_CB_local_prev::Array{Float64, 2}
    P_BESS_local_prev::Array{Float64, 2}
end

function setupAsInPGM(AgentCB::Array{Any,1}, N::Int64, T::Int64)

    # initialize
    async_pgm_u = Array{Any}(N)
    async_pgm_y = Array{Any}(N)
    async_pgm_x = Array{Any}(N)
    async_pgm_v = Array{Any}(N)

    for i = 1:N
        async_pgm_u[i] = zeros(AgentCB[i].system.Nu,T)
        async_pgm_y[i] = zeros(AgentCB[i].system.Ny,T)
        async_pgm_x[i] = zeros(AgentCB[i].system.Nx,T)
        async_pgm_v[i] = zeros(T)  # agent's consumption
    end
    async_pgm_P_CB = zeros(T,N)  # agents' consumptions
    async_pgm_P_CB_prev = ones(T,N)
    async_pgm_P_BESS = zeros(T)  # battery's consumption
    async_pgm_P_BESS_prev = ones(T)
    async_pgm_RES = [0.0]

    # matrices containing the last read of the global vector for each agent
    async_pgm_P_CB_local = zeros(N*T,N+1)
    async_pgm_P_BESS_local = zeros(T,N+1)
    async_pgm_P_CB_local_prev = zeros(N*T,N+1)
    async_pgm_P_BESS_local_prev = zeros(T,N+1)

    # pgm = PGM(pgm_u, pgm_y, pgm_x, pgm_v, pgm_P_CB, pgm_P_CB_prev, pgm_P_BESS, pgm_P_BESS_prev)

    return async_pgm_u, async_pgm_y, async_pgm_x, async_pgm_v, async_pgm_P_CB, async_pgm_P_CB_prev, async_pgm_P_BESS, async_pgm_P_BESS_prev, async_pgm_RES, async_pgm_P_CB_local, async_pgm_P_BESS_local, async_pgm_P_CB_local_prev, async_pgm_P_BESS_local_prev
end

function AlgoAsInPGM(coordinate::Bool, varying_β::Bool, as_in_pgm::AsInPGM, N::Int64, T::Int64, Tend::Float64,
  UpdateQueue::Dict{Int64,Float64},
  ExecTimesMean::Array{Float64,1}, ExecTimesVar::Array{Float64,1},
  optimizationModelCB::Array{Any,1}, optimizationModelBESS::ProxBESS,
  agentBESS::BESS, agentCB::Array{Any,1},
  master::Master, p_BESS_opt::Array{Float64,1}, v_p_CB_opt::Array{Float64,1}, η::Float64, β::Float64, δ::Float64, α::Float64)

    tempCB   = zeros(T,N)
    tempBESS = zeros(T)

    # intervals = Int(Tend/0.2 + 2)
    NoUpdates = zeros(N+1)
    LastUpdateIndex = ones(Int64,N+1)
    k = 1
    T2 = [0.0]

    count = 1

    # main loop
    while T2[end] < Tend

        # choose the agent that's next in line to update
        minIdx = indmin(collect(values(UpdateQueue)))
        TimeNextUpdate = collect(values(UpdateQueue))[minIdx]
        ik = collect(keys(UpdateQueue))[minIdx]
            # display("Update queue: $(collect(values(UpdateQueue)))")
            # display("Time of the next update: $(TimeNextUpdate)")
            # display("Next index to update: $(ik)")
        NoUpdates[ik] += 1

        if ik == N+1

            # Master computes residual
            res = computeResidual(master, α, as_in_pgm.P_BESS_local[:,N+1], as_in_pgm.P_CB_local[:,N+1])

            # BESS proximal step
            as_in_pgm.P_BESS_prev = as_in_pgm.P_BESS
            point = computeProxPointBESS(agentBESS, δ, β, as_in_pgm.P_BESS_local[:,N+1], as_in_pgm.P_BESS_local_prev[:,N+1], res)
            # tempBESS = getvalue(solveProxBESS(agentBESS, T, point, SOC0, δ))
            tempBESS = getvalue(solveProxBESS(optimizationModelBESS, point, sparse(agentBESS.data.Γ)))
            # tempBESS = evaluate(solveProxBESS_Convex(agentBESS, T, point, P0_half, q0, SOC0, δ))
            as_in_pgm.P_BESS += η * (tempBESS - as_in_pgm.P_BESS)

            if !coordinate
                # update the building agents with the lastly computed prox
                for j in 1:N
                    as_in_pgm.P_CB_prev[:,j] = as_in_pgm.P_CB[:,j]
                    as_in_pgm.P_CB[:,j] = (1-η) * as_in_pgm.P_CB[:,j] + η *  tempCB[:,j]
                end
            end

            # update battery's knowledge (to be used in the next calculation)
            as_in_pgm.P_CB_local[:,N+1] = vec(as_in_pgm.P_CB)
            as_in_pgm.P_BESS_local_prev[:,N+1] = as_in_pgm.P_BESS_prev
            as_in_pgm.P_BESS_local[:,N+1] = as_in_pgm.P_BESS

        else

            # Master computes residual
            res = computeResidual(master, α, as_in_pgm.P_BESS_local[:,ik], as_in_pgm.P_CB_local[:,ik])

            # building proximal step
            as_in_pgm.P_CB_prev[:,ik] = as_in_pgm.P_CB[:,ik]
            point = computeProxPointCB(agentCB[ik], δ, β, as_in_pgm.P_CB_local[(ik-1)*T+1:ik*T,ik], as_in_pgm.P_CB_local_prev[(ik-1)*T+1:ik*T,ik], res)
            # as_in_pgm.u[ik], as_in_pgm.y[ik], as_in_pgm.x[ik], as_in_pgm.v[ik] = solveProxCB(agentCB[ik], Tref, T, X0[ik], W[ik], point, night_times, day_times, δ)
            as_in_pgm.u[ik], as_in_pgm.y[ik], as_in_pgm.x[ik], as_in_pgm.v[ik] = solveProxCB(optimizationModelCB[ik], point, sparse(agentCB[ik].data.Γ))
            tempCB[:,ik] = getvalue(as_in_pgm.v[ik])

            # as_in_pgm.u[ik], as_in_pgm.y[ik], as_in_pgm.x[ik], as_in_pgm.v[ik] = solveProxCB_Convex(agentCB[ik], Tref, T, X0[ik], W[ik], point, night_times, day_times, δ)
            # tempCB[:,ik] = evaluate(as_in_pgm.v[ik])

            as_in_pgm.P_CB[:,ik] += η * (tempCB[:,ik] - as_in_pgm.P_CB[:,ik])

            if !coordinate
                # update the rest of the components
                for j in 1:N+1
                    if j == N+1
                      as_in_pgm.P_BESS_prev = as_in_pgm.P_BESS
                      as_in_pgm.P_BESS += η * (tempBESS - as_in_pgm.P_BESS)
                    elseif ((j!=(N+1)) && (j!=ik))
                      as_in_pgm.P_CB_prev[:,j] = as_in_pgm.P_CB[:,j]
                      as_in_pgm.P_CB[:,j] += η * (tempCB[:,j] - as_in_pgm.P_CB[:,j])
                    end
                end
            end

            # update building's knowledge
            as_in_pgm.P_BESS_local[:,ik] = as_in_pgm.P_BESS
            as_in_pgm.P_CB_local_prev[:,ik] = vec(as_in_pgm.P_CB_prev)
            as_in_pgm.P_CB_local[:,ik]   = vec(as_in_pgm.P_CB)

        end

        # residual
        push!(as_in_pgm.RES, norm(as_in_pgm.P_BESS - p_BESS_opt)/norm(p_BESS_opt) + norm(vec(as_in_pgm.P_CB)-v_p_CB_opt)/norm(v_p_CB_opt))

        # save time values
        push!(T2, TimeNextUpdate)
        if mod(k, 100) == 0
            display("Update queue: $(collect(values(UpdateQueue)))")
            display("Time of the next update: $(TimeNextUpdate)")
            display("Next index to update: $(ik)")
            display("Time is: $(T2[end])")
            display("Residual: $(as_in_pgm.RES[end])")
            count += 1
            if varying_β 
                β = 1 - 2/(count+1)
                display("New β is: $(β)")
            end
        end

        # update the queue
        UpdateQueue[ik] = TimeNextUpdate + ExecTimesMean[ik] + ExecTimesVar[ik]*randn()
        k += 1

    end

    return as_in_pgm, T2, NoUpdates

end

end
