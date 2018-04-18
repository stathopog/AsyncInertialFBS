# AsyncInertialFBS
Julia code that implements an asynchronous version of the Proximal Gradient Method (PGM) (instant of the Forward-Backward Splitting method) and compares with other implementations of the algorithm by solving a load sharing microgrid problem

## What it does
Julia notebook that solves a load-sharing problem in a microgrid setting comprising controllable buildings (CBs) and a battery electical storage system (BESS).
The goal is to collaboratively track a generated signal while complying to local constraints.
Each agent (CB or BESS) solves a proximal minimization problem while a central coordinator collects the solutions and transmits new directions back to the agents. 
Due to the heterogeneity of the agents' computational times, an asynchronous version of PGM is employed and compared to other implementations.  

## How to use it
The main file is the  jupyter notebook. Imports the following files:
* **battery.jl:**              the BESS module (constructor)
* **building.jl:**             the CB module (constructor for three building types)
* **master.jl:**               the coordinator module (constructor)
* **parseData.jl:**            the original building models have been generated in MATLAB using the OpenBuild package; module imports them from the *data* folder
* **agentCB.jl:**              aggregator of buidlings that setups the local proximal minimization problems
* **agentBESS.jl:**            battery (BESS) module
* **pgm.jl:**                        PGM module
* **async_inert_pgm.jl:**   asynchronous and inertial PGM module

All optimization problems are solved using JuMP with the Gurobi solver.

## References
1. G. Stathopoulos and C. N. Jones, "An Inertial Parallel and Asynchronous Forward-Backward Iteration for Distributed Convex Optimization", arXiv:1706.00088
