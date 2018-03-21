# Battery agent data

#=
Sets up the battery node.
The module includes the following components:

Optimization:                        type containing the several parameters that
                                     enter the optimization problem, e.g., hexogenous signals
                                     (disturbances/predictions/references),
                                     algorithmic constants etc.

BESS:                                type describing the battery entity; dynamics and optimization parameters

ProxBESS:                            type describing the optimization model (proximal minimization problem) that is going to
                                     be repeatedly solved later. Constitutes of JuMP objects

setupProxCB:                         function tha parses the proximal minimization problem

computeProxPointBESS:                function that computes the point at which the proximal minimization will be solved

solveProxBESS:                       function that solves the proximal minimization problem
                                     associated to the battery


@author : giorgos stathopoulos

@date : 2017-06-20
=#

module AgentsBESS

include("battery.jl")
using Batteries, Gurobi, JuMP

export BESS, computeProxPointBESS, ProxBESS, solveProxBESS

type Optimization
    Γ::Array{Float64,2}        # preconditioner
    c_bess::Array{Float64, 1}  # affine objective part
end

type BESS
    system::Battery
    data::Optimization
end

type ProxBESS
    batteryModel::JuMP.Model
    p_BESS::Array{JuMP.Variable,1}
    obj::JuMP.GenericQuadExpr{Float64,JuMP.Variable}
end

function setupBESS(battery::Battery, SOC0::Float64, SOCref::Float64, Γ::Array{Real,2} )

    c = battery.AA_bess*SOC0 - SOCref
    data = Optimization(Γ, c)

    system = Battery(battery.A, battery.Bu, battery.Bw, battery.AA_bess, battery.BB_bess,
                     battery.nx, battery.nu, battery.pMin, battery.pMax, battery.SoCMin, battery.SoCMax)

    return system, data
end

function computeProxPointBESS(agent::BESS, δ::Float64, β::Float64, p::Array{Float64, 1}, p_prev::Array{Float64, 1}, res::Array{Float64, 1})

    # compute point where proximal opearator is evaluated
    v = p - agent.data.Γ\res + β*(p-p_prev)

    return v
end

function setupProxBESS(agent::BESS,
                       T::Int64,
                       # point::Array{Float64, 1},
                       SOC0::Float64, δ::Float64)

    # initialize model
    batteryModel = Model(solver=GurobiSolver(OutputFlag=0))

    @variable(batteryModel, p_BESS[1:T]) # total power consumption

    # inequality constraints
    for t = 1:T
       @constraint(batteryModel, p_BESS[t] >= agent.system.pMin)
       @constraint(batteryModel, p_BESS[t] <= agent.system.pMax)
       @constraint(batteryModel, agent.system.AA_bess[t]*SOC0 + dot(agent.system.BB_bess[t,:], p_BESS) <= agent.system.SoCMax)
       @constraint(batteryModel, agent.system.AA_bess[t]*SOC0 + dot(agent.system.BB_bess[t,:], p_BESS) >= agent.system.SoCMin)
    end

    obj = dot(p_BESS, 0.5*δ*(agent.system.BB_bess'*agent.system.BB_bess)*p_BESS) + δ*dot(p_BESS, (agent.system.BB_bess'*agent.data.c_bess)) + 0.5*δ*dot(agent.data.c_bess, agent.data.c_bess)

    @objective( batteryModel, Min, obj )

    solve(batteryModel)

    return ProxBESS(batteryModel, p_BESS, obj)
end


function solveProxBESS(optimizationModel::ProxBESS,
                       point::Array{Float64,1},
                       H::SparseMatrixCSC{Float64,Int64})
    obj_prox = 0.5*vecdot(optimizationModel.p_BESS-point, H*(optimizationModel.p_BESS-point))
    @objective( optimizationModel.batteryModel, Min, optimizationModel.obj + obj_prox )

    solve(optimizationModel.batteryModel)

    return optimizationModel.p_BESS
end

# function solveProxBESS_Convex(agent::BESS, T::Int64, point::Array{Float64, 1}, P0_half::Array{Float64, 2}, q0::Array{Float64, 1}, SOC0::Float64, δ::Float64)

#     # initialize model
#     solver = GurobiSolver(OutputFlag=0)

#     p_BESS = Variable(T) # total power consumption
#     OBJ_Track  = Variable(1) # epigraph form for SOC tracking objective
#     OBJ_Prox  = Variable(1) # epigraph form for proximal regularization term

#     problem = minimize( OBJ_Track + sum(OBJ_Prox) )

#     # inequality constraints - BESS
#     problem.constraints += (p_BESS <= agent.system.pMax)
#     problem.constraints += (p_BESS >= agent.system.pMin)
#     problem.constraints += (agent.system.AA_bess*SOC0 + agent.system.BB_bess*p_BESS <= agent.system.SoCMax)
#     problem.constraints += (agent.system.AA_bess*SOC0 + agent.system.BB_bess*p_BESS >= agent.system.SoCMin)

#     # SOC constraints for objective
#     problem.constraints += ( norm(P0_half*p_BESS + P0_half\q0) <= OBJ_Track )

#     # proximal term
#     problem.constraints += ( norm(sqrt(1/2)*sqrtm(agent.data.Γ)\p_CB + sqrt(2)*sqrtm(agent.data.Γ)*(-point) ) <= OBJ_Prox )

#     solve!(problem, solver)

#     return p_BESS
# end




end
