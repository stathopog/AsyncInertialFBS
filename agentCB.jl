# Building agent data

#=
Sets up an aggregator comprising buildings.
The module includes the following components:

Consumption:        type containing four consumption-related characteristics of the buildings

StateSpace:         type w/ the state-space description of the aggregation

Constraints:        type w/ the constraints (assumed to be the same for
                    all types of buildings)

Optimization:       type w/ the optimization parameters Γ (stepsizes matrix) and α (cost regularizer)

Aggregator:         type; the actual agent that constitutes of the types defined above

ProxCB:             type describing the optimization model (proximal minimization problem) that is going to
                    be repeatedly solved later. Constitutes of JuMP objects

setupProxCB:        function tha parses the proximal minimization problem

solveProx:          function that solves the proximal minimization problem
                    associated to the aggregator

computeProxPointCB: function that computes the point at which the proximal minimization will be solved

@author : giorgos stathopoulos

@date : 2017-05-15
=#

module AgentsCB

include("building.jl")
using Buildings, Gurobi, JuMP

export Consumption, StateSpace, Constraints, Aggregator, computeProxPointCB, ProxCB, solveProxCB

type Consumption
    Uss::Array{Float64, 2}                   # nominal input to preserve baseline consumption
    baseline_thermal::Array{Float64, 2}      # thermal baseline consumption
    COP::Float64                             # Coefficient Of Performance (electrical -> thermal)
    baseline_electrical::Array{Float64, 2}   # electrical baseline consumption
end

type StateSpace
    A::Array{Float64, 2}
    Bu::Array{Float64, 2}
    Bw::Array{Float64, 2}
    C::Array{Float64, 2}
    Nx::Int64
    Nu::Int64
    Nw::Int64
    Ny::Int64
end

type Constraints
    inputMin::Matrix
    inputMax::Matrix
    outputDayMin::Matrix
    outputDayMax::Matrix
    outputNightMin::Matrix
    outputNightMax::Matrix
end

type Optimization
    Γ::Array{Float64, 2}  # preconditioner
    α::Float64  # penalty to dispatch plan
end

type ProxCB
    agentModel::JuMP.Model
    U::Array{JuMP.Variable,2}
    X::Array{JuMP.Variable,2}
    Y::Array{JuMP.Variable,2}
    p_CB::Array{JuMP.Variable,1}
    obj::JuMP.GenericQuadExpr{Float64,JuMP.Variable}
end

type Aggregator
    subsystems::Int16
    consumption::Consumption
    system::StateSpace
    constraints::Constraints
    data::Optimization
end

function setupAggregator(No::Int64, building::Building, Γ::Array{Float64, 2}, α::Float64)

    A    = Array{Float64}(0,0)
    Bu   = Array{Float64}(0,0)
    Bw   = Array{Float64}(0,0)
    C    = Array{Float64}(0,0)

    # consumptions
    Uss                 = kron(ones(No), building.Uss)
    baseline_thermal    = repmat(building.baseline, 1, No)
    COP                 = building.COP
    baseline_electrical = baseline_thermal' / (1000*COP)

    consume = Consumption(Uss, baseline_thermal, COP, baseline_electrical)

    # dynamics
    A  = kron(eye(No),building.A)
    Bu = kron(eye(No),building.Bu)
    Bw = kron(eye(No),building.Bw)
    C  = kron(eye(No),building.C)

    # dimensions
    Nx = size(A,1)
    Nu = size(Bu,2)
    Nw = size(Bw,2)
    Ny = size(C,1)

    system = StateSpace(A, Bu, Bw, C, Nx, Nu, Nw, Ny)

    # constraints
    inputMin = building.inputMin
    inputMax = building.inputMax
    outputDayMin = building.outputDayMin
    outputDayMax = building.outputDayMax
    outputNightMin = building.outputNightMin
    outputNightMax = building.outputNightMax
    constraint = Constraints(inputMin, inputMax, outputDayMin, outputDayMax, outputNightMin, outputNightMax)

    # optimization data
    Γ  = kron(ones(No), Γ)
    data = Optimization(Γ, α)

    return consume, system, constraint, data
end


function computeProxPointCB(agent::Aggregator, δ::Float64, β::Float64, p::Array{Float64, 1}, p_prev::Array{Float64, 1}, res::Array{Float64, 1})

    # compute point where proximal opearator is evaluated
    v = p - agent.data.Γ\(res + δ*(p-agent.consumption.baseline_electrical)) + β*(p-p_prev)

    return v
end


function setupProxCB(agent::Aggregator, Tref::Float64, T::Int64,
                 x0::Array{Float64, 1}, w::Array{Float64, 2},
                 night_times::Array{Int64, 1}, day_times::Array{Int64, 1}, δ::Float64)

    # initialize model
    agentModel = Model(solver=GurobiSolver(OutputFlag=0))

    @variable(agentModel, U[1:agent.system.Nu,1:T])
    @variable(agentModel, X[1:agent.system.Nx,1:T])
    @variable(agentModel, Y[1:agent.system.Ny,1:T])
    @variable(agentModel, p_CB[1:T]) # total power consumption

    # dynamics & equality contraints
    for i = 1:agent.system.Nx
        @constraint(agentModel, (X[i,1] == vecdot(agent.system.A[i,:], x0) + vecdot(agent.system.Bu[i,:], U[:,1]) +
                                          vecdot(agent.system.Bw[i,:], w[:,1]) ) )
    end
    for j = 1:agent.system.Ny
        @constraint(agentModel, Y[j,1] == vecdot(agent.system.C[j,:], X[:,1]) )
    end

    obj = 0.0
    for t = 1:T-1
        for i = 1:agent.system.Nx
            @constraint(agentModel, (X[i,t+1] == vecdot(agent.system.A[i,:], X[:,t]) + vecdot(agent.system.Bu[i,:], U[:,t+1]) +
                                                vecdot(agent.system.Bw[i,:], w[:,t+1]) ) )
        end

        for j = 1:agent.system.Ny
            @constraint(agentModel, Y[j,t+1] == vecdot(agent.system.C[j,:], X[:,t+1]) )
        end

        @constraint(agentModel, p_CB[t] == sum(U[:,t])/(1000*agent.consumption.COP))
    end
    @constraint(agentModel, p_CB[T] == sum(U[:,T])/(1000*agent.consumption.COP))

    # inequality constraints
    for t = 1:T
        if ((t-1) in night_times)
            for j = 1:agent.system.Ny
               @constraint(agentModel, Y[j,t]  >= agent.constraints.outputNightMin[j])
               @constraint(agentModel, Y[j,t]  <= agent.constraints.outputNightMax[j])
            end
        elseif ((t-1) in day_times)
            for j = 1:agent.system.Ny
               @constraint(agentModel, Y[j,t]  >= agent.constraints.outputDayMin[j])
               @constraint(agentModel, Y[j,t]  <= agent.constraints.outputDayMax[j])
            end
        end

        for i = 1:agent.system.Nu # for all rows do the following
            @constraint(agentModel, 0*agent.constraints.inputMin[i] <= U[i,t])
            @constraint(agentModel, U[i,t] <= agent.constraints.inputMax[i])
        end
    end

    # reference temperature tracking
    for t = 1:T
        obj = obj + 0.5*δ*vecdot(Y[:,t] - Tref, Y[:,t] - Tref)
    end

    @objective( agentModel, Min, obj )

    solve(agentModel)

    return ProxCB(agentModel, U, X, Y, p_CB, obj)
end


function solveProxCB(optimizationModel::ProxCB,
                     point::Array{Float64,2},
                     H::SparseMatrixCSC{Float64,Int64})

    obj_prox = 0.5*vecdot(optimizationModel.p_CB-point, H*(optimizationModel.p_CB-point))
    @objective( optimizationModel.agentModel, Min, optimizationModel.obj + obj_prox )

    solve(optimizationModel.agentModel)

    return optimizationModel.U, optimizationModel.Y, optimizationModel.X, optimizationModel.p_CB
end


# function solveProxCB_Convex(agent::Aggregator, Tref::Float64, T::Int64,
#                  x0::Array{Float64, 1}, w::Array{Float64, 2},
#                  point::Array{Float64, 2},
#                  night_times::Array{Int64, 1}, day_times::Array{Int64, 1}, δ::Float64)

#     # initialize model
#     solver = GurobiSolver(OutputFlag=0)

#     U    = Variable(agent.system.Nu,T)
#     X    = Variable(agent.system.Nx,T)
#     Y    = Variable(agent.system.Ny,T)
#     p_CB = Variable(T) # total power consumption
#     OBJ_Temp  = Variable(1) # epigraph form for temperature tracking objective
#     OBJ_Prox  = Variable(1) # epigraph form for proximal regularization term

#     problem = minimize( OBJ_Temp + OBJ_Prox )

#     # dynamics & equality contraints
#     problem.constraints += ( X[:,1] == agent.system.A*x0 + agent.system.Bu*U[:,1] + agent.system.Bw*w[:,1] )
#     problem.constraints += ( Y[:,1] == agent.system.C*X[:,1] )

#     problem.constraints += ( X[:,2:T] == agent.system.A*X[:,1:T-1] + agent.system.Bu*U[:,2:T] + agent.system.Bw*w[:,2:T] )
#     problem.constraints += ( Y[:,2:T] == agent.system.C*X[:,2:T] )

#     temp = sum(U,1)
#     for t in 1:T
#         problem.constraints += ( p_CB[t] == temp[t]/(1000*agent.consumption.COP) )
#     end

#     # inequality constraints
#     [problem.constraints += Y[:,t+1] >= agent.constraints.outputNightMin for t in intersect(collect(1:T)-1,night_times) ]
#     [problem.constraints += Y[:,t+1] <= agent.constraints.outputNightMax for t in intersect(collect(1:T)-1,night_times) ]
#     [problem.constraints += Y[:,t+1] >= agent.constraints.outputDayMin for t in intersect(collect(1:T)-1,day_times) ]
#     [problem.constraints += Y[:,t+1] <= agent.constraints.outputDayMax for t in intersect(collect(1:T)-1,day_times) ]
#     [problem.constraints += 0*agent.constraints.inputMin <= U[:,t] for t = 1:T]
#     [problem.constraints += agent.constraints.inputMax >= U[:,t] for t = 1:T]

#     # reference temperature tracking
#     problem.constraints += ( norm(sqrt(δ/2)*vec(Y) + sqrt(2/δ)*(-δ*0.5*Tref*ones(agent.system.Ny*T,1)) ) <= OBJ_Temp )

#     # proximal term
#     problem.constraints += ( norm(sqrt(1/2)*(1/sqrt(agent.data.Γ[1,1]))*p_CB + sqrt(2)*sqrt(agent.data.Γ[1,1])*(-point) ) <= OBJ_Prox )

#     solve!(problem, solver)

#     return U, Y, X, p_CB
# end

end
