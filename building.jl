# Building model definition

#=
Constructor of a linear building model
for MPC of the form

    x(t+1) = A*x(t) + Bu*u(t) + Bw*w(t)

    y(t) = C*x(t)

with box constraints on u(t) for the
consumption
and
box contraints on y(t) preventing deviations
from the desired temperature.

The module includes the following components:

Building:          type containing the full description of a building (state-space, dimensions, constraints)

parseBuildings:   function that reads in MATLAB data for the three building types (small, medium, large)

@author : giorgos stathopoulos

@author : giorgos stathopoulos

@date : 2017-05-14
=#

module Buildings

using DataFrames

export Building

type Building
    class::AbstractString
    Uss::Matrix
    baseline::Matrix
    A::Matrix
    Bu::Matrix
    Bw::Matrix
    C::Matrix
    inputMin::Matrix
    inputMax::Matrix
    outputDayMin::Matrix
    outputDayMax::Matrix
    outputNightMin::Matrix
    outputNightMax::Matrix
    COP::Float64
end

function parseBuildings(build_instance::Array{Building,1})

    # Read in buildings
    for i in 1:3
        U_ss     = convert(Array,readtable(string("data/U_ss_", "$(i)", ".dat"), separator = ',', header = false))
        baseline = collect(convert(Array,readtable(string("data/base_", "$(i)", ".dat"), separator = ',', header = false)))

        # read in dynamics
        A  = convert(Array,readtable(string("data/A_", "$(i)", ".dat"), separator = ',', header = false))
        Bu = convert(Array,readtable(string("data/Bu_", "$(i)", ".dat"), separator = ',', header = false))
        Bw = convert(Array,readtable(string("data/Bw_", "$(i)", ".dat"), separator = ',', header = false))
        C  = convert(Array,readtable(string("data/C_", "$(i)", ".dat"), separator = ',', header = false))

        # read in constraints
        inputMin  = collect(convert(Array,readtable(string("data/inputMin_", "$(i)", ".dat"), separator = ',', header = false)))
        inputMax  = collect(convert(Array,readtable(string("data/inputMax_", "$(i)", ".dat"), separator = ',', header = false)))
        outputDayMin = collect(convert(Array,readtable(string("data/outputDayMin_", "$(i)", ".dat"), separator = ',', header = false)))
        outputDayMax = collect(convert(Array,readtable(string("data/outputDayMax_", "$(i)", ".dat"), separator = ',', header = false)))
        outputNightMin = collect(convert(Array,readtable(string("data/outputNightMin_", "$(i)", ".dat"), separator = ',', header = false)))
        outputNightMax = collect(convert(Array,readtable(string("data/outputNightMax_", "$(i)", ".dat"), separator = ',', header = false)))
        if i == 1
            class = "small"
        elseif i == 2
            class = "medium"
        elseif i == 3
            class = "large"
        end

        COP = 3.0

        build_instance[i] = Building(class, U_ss, baseline, A, Bu, Bw, C, inputMin, inputMax, outputDayMin, outputDayMax, outputNightMin, outputNightMax, COP)
        println("\n$(size(build_instance[i].Bu))\n")
    end
    return build_instance
end

end
