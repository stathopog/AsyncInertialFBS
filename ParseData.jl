# ParseData

#=
Parse buildings and other data

parseSignals:          reads in initial states and disturbance signal

parsePreconditioners:  reads in the preconditioning matrices (stepsizes)

@author : giorgos stathopoulos

@date : 2017-06-19
=#

module ParseData

# include("building.jl")
# include("master.jl")
using DataFrames

function parseSignals(x0::Array{Any}, w::Array{Any}, ref::Vector)

    ref = collect(convert(Array,readtable(string("data/ref", ".dat"), separator = ',', header = false)))
    ref = squeeze(ref,2)

    for i in 1:3
        # read in initial state and noise vector
        x0[i] = collect(convert(Array,readtable(string("data/x0_", "$(i)", ".dat"), separator = ',', header = false)))
        x0[i] = squeeze(x0[i],2)
        w[i] = collect(convert(Array,readtable(string("data/w_", "$(i)", ".dat"), separator = ',', header = false)))
    end
    return x0, w, ref
end

function parsePreconditioners(Q_BESS::Array{Any}, Q_CB::Array{Any})

    Q_BESS  = convert(Array,readtable(string("data/Gamma_1", ".dat"), separator = ',', header = false))
    Q_CB    = convert(Array,readtable(string("data/Gamma_2", ".dat"), separator = ',', header = false))
    Q_BESS  = convert(Array{Real,2}, Q_BESS)
    Q_CB    = convert(Array{Float64,2}, Q_CB)

    return Q_BESS, Q_CB
end


# function parseBattery(battery_instance::Array{Battery,1}, N::Int64)

#     # Read in battery

#     # read in dynamics
#     A  = convert(Array,readtable(string("data/A_master", ".dat"), header = false))
#     Bu = convert(Array,readtable(string("data/Bu_master", ".dat"), separator = ',', header = false))
#     Bw = convert(Array,readtable(string("data/Bw_master", ".dat"), header = false))

#     # read in dynamics in condensed form
#     A_batt  = convert(Array,readtable(string("data/Abatt_master", ".dat"), separator = ',', header = false))
#     B_batt  = convert(Array,readtable(string("data/Bbatt_master", ".dat"), separator = ',', header = false))
#     # MSoC0 = convert(Array,readtable(string("data/Mx0_master", ".dat"), separator = ',', header = false))
#     # Ms = convert(Array,readtable(string("data/Ms_master", ".dat"), separator = ',', header = false))
#     # CC = convert(Array,readtable(string("data/CC_master", ".dat"), separator = ',', header = false))

#     # read in dimensions
#     nx = 1
#     nu = 1
#     # nw = 1

#     # # read in constraints
#     pMin  = -5.0
#     pMax  = 5.0
#     SoCMin  = 0.0
#     SoCMax  = 40.0

#     battery_instance = Battery(A, Bu, Bw, A_batt, B_batt, nx, nu, pMin, pMax, SoCMin, SoCMax)

#     return battery_instance
# end


# function parseBuildings(build_instance::Array{Building,1})

#     # Read in buildings
#     for i in 1:3
#         U_ss     = convert(Array,readtable(string("data/U_ss_", "$(i)", ".dat"), separator = ',', header = false))
#         baseline = collect(convert(Array,readtable(string("data/base_", "$(i)", ".dat"), separator = ',', header = false)))

#         # read in dynamics
#         A  = convert(Array,readtable(string("data/A_", "$(i)", ".dat"), separator = ',', header = false))
#         Bu = convert(Array,readtable(string("data/Bu_", "$(i)", ".dat"), separator = ',', header = false))
#         Bw = convert(Array,readtable(string("data/Bw_", "$(i)", ".dat"), separator = ',', header = false))
#         C  = convert(Array,readtable(string("data/C_", "$(i)", ".dat"), separator = ',', header = false))

#         # read in dynamics in condensed form
#         Mu  = convert(Array,readtable(string("data/Mu_", "$(i)", ".dat"), separator = ',', header = false))
#         Mx0 = convert(Array,readtable(string("data/Mx0_", "$(i)", ".dat"), separator = ',', header = false))
#         Mw = convert(Array,readtable(string("data/Mw_", "$(i)", ".dat"), separator = ',', header = false))
#         out_Mu  = convert(Array,readtable(string("data/out_Mu_", "$(i)", ".dat"), separator = ',', header = false))
#         out_Mx0 = convert(Array,readtable(string("data/out_Mx0_", "$(i)", ".dat"), separator = ',', header = false))
#         out_Mw = convert(Array,readtable(string("data/out_Mw_", "$(i)", ".dat"), separator = ',', header = false))

#         # read in constraints
#         inputMin  = collect(convert(Array,readtable(string("data/inputMin_", "$(i)", ".dat"), separator = ',', header = false)))
#         inputMax  = collect(convert(Array,readtable(string("data/inputMax_", "$(i)", ".dat"), separator = ',', header = false)))
#         outputDayMin = collect(convert(Array,readtable(string("data/outputDayMin_", "$(i)", ".dat"), separator = ',', header = false)))
#         outputDayMax = collect(convert(Array,readtable(string("data/outputDayMax_", "$(i)", ".dat"), separator = ',', header = false)))
#         outputNightMin = collect(convert(Array,readtable(string("data/outputNightMin_", "$(i)", ".dat"), separator = ',', header = false)))
#         outputNightMax = collect(convert(Array,readtable(string("data/outputNightMax_", "$(i)", ".dat"), separator = ',', header = false)))
#         if i == 1
#             class = "small"
#         elseif i == 2
#             class = "medium"
#         elseif i == 3
#             class = "large"
#         end
#         build_instance[i] = Building(class, U_ss, baseline, A, Bu, Bw, C, Mu, Mx0, Mw, out_Mu, out_Mx0, out_Mw,
#          inputMin, inputMax, outputDayMin, outputDayMax, outputNightMin, outputNightMax)
#         println("\n$(size(build_instance[i].Bu))\n")
#     end
#     return build_instance
# end


end

