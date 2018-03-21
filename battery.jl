# Battery model definition

#=
Constructor of a linear battery model
for MPC of the form

    SoC(t+1) = a*SoC(t) + p_batt(t) + s(ϵ_unc(t),t)

with box constraints on p_batt(t).

The model can be rewritten as

    SoC(t+1) = A_batt*SoC(0) + B_batt*p_CB(t) + B_ϵ,batt*s(ϵ_unc(t),t)

    D*p_CB < d,

using the equation

    p_batt(t) = -sum_i (\hat{p}_i,CB(t) + p_i,CB(t) + ϵ_unc(t).

The module includes the following components:

Battery:        type containing the full description of the battery (state-space, dimensions, constraints)

parseBattery:   function that reads in MATLAB data that generated the model

condense:       function that condenses the model (inputs-only description)

@author : giorgos stathopoulos

@date : 2017-06-20
=#

module Batteries

using DataFrames

export Battery

type Battery
    A::Array{Float64, 2}
    Bu::Array{Float64, 2}
    Bw::Array{Float64, 2}
    AA_bess::Array{Float64, 1}
    BB_bess::Array{Float64, 2}
    nx::Int64
    nu::Int64
    pMin::Float64
    pMax::Float64
    SoCMin::Float64
    SoCMax::Float64
end

function parseBattery(battery_instance::Array{Battery,1}, N::Int64, T::Int64)

    # Read in battery

    # read in dynamics
    A  = convert(Array,readtable(string("data/A_bess", ".dat"), header = false))
    Bu = convert(Array,readtable(string("data/Bu_bess", ".dat"), separator = ',', header = false))
    Bw = convert(Array,readtable(string("data/Bw_bess", ".dat"), header = false))

    # read in dynamics in condensed form
    AA_bess  = convert(Array,readtable(string("data/AA_bess", ".dat"), separator = ',', header = false))
    AA_bess = squeeze(AA_bess,2)
    BB_bess  = convert(Array,readtable(string("data/BB_bess", ".dat"), separator = ',', header = false))

    # read in dimensions
    nx = 1
    nu = 1
    # nw = 1

    # read in constraints
    SoCMin  = 0.0
    SoCMax  = 1.0
    pMin    = -SoCMax / 5
    pMax    = SoCMax / 5

    battery_instance = Battery(A, Bu, Bw, AA_bess, BB_bess, nx, nu, pMin, pMax, SoCMin, SoCMax)

    return battery_instance
end

function condense(battery_instance::Battery, T::Int64)
    A_bess = battery_instance.A
    for i = 2:T
        A_bess = [A_bess; battery_instance.A^i]
    end
    bar_B = battery_instance.Bu
    for i = 1:T-1
        bar_B = [bar_B; battery_instance.A^i*battery_instance.Bu]
    end
    B_bess = [bar_B zeros(T*battery_instance.Nx,(T-1)*battery_instance.Nu)]
    for i = 1:T-1
        B_bess[:,i*battery_instance.Nu+1:(i+1)*battery_instance.Nu] = [ zeros(i*battery_instance.Nx,battery_instance.Nu); bar_B[1:end-i*battery_instance.Nx,:] ]
    end

    # returns the A and B battery matrices
    return A_bess, B_bess
end

end
