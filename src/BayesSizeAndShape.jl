module BayesSizeAndShape
 
##### Packages
using Distributions, Random
using LinearAlgebra, PDMats
using ProgressMeter
using Kronecker
using DataFrames
using CategoricalArrays
using StatsModels
using ToggleableAsserts

##### Include
include(joinpath("types.jl"))
include(joinpath("general_functions.jl"))
include(joinpath("mcmc.jl"))
include(joinpath("sampler.jl"))
include(joinpath("wrappers.jl"))
include(joinpath("deprecated.jl"))

##### Functions
 export
    NoPriorBeta,
    NoPriorSigma,
    compute_ss_from_pre,
    KeepReflection,
    DoKeepReflection,
    GeneralSigma,
    standardize_reg,
    SizeAndShapeMCMC,
    Valuep2,
    compute_designmatrix,
    compute_dessignmatrix,
    compute_helmertized_configuration

end 
