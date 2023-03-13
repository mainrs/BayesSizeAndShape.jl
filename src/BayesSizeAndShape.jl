module BayesSizeAndShape
 

#using Pkg
#Pkg.activate("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/gitrepo/BayesSizeAndShape")

##### Packages
# using MKL
using Distributions, Random
using LinearAlgebra, PDMats
using ProgressMeter
using Kronecker
using DataFrames
using StatsModels
using ToggleableAsserts
# using ThreadsX, SparseArrays
# using ThreadedSparseArrays
# using SparseMatricesCSR
# using ThreadedSparseCSR
# # using Impute
# # using SpecialFunctions

##### Include
include(joinpath("types.jl"))
include(joinpath("general_functions.jl"))
include(joinpath("mcmc.jl"))
include(joinpath("sampler.jl"))
include(joinpath("wrappers.jl"))
# include(joinpath("func_covariance_new.jl"))
# include(joinpath("func_project.jl"))
# include(joinpath("new_func_project.jl"))
# include(joinpath("func_mcmc.jl"))


# include(joinpath("func_adapt.jl"))
# include(joinpath("sampler.jl"))
# include(joinpath("sampler_new.jl"))

# # include(joinpath("sampler_notraset.jl")) 
# # include(joinpath("func_adapt_notraset.jl")) 


# include(joinpath("predictw.jl"))



##### Functions
 export
    NoPriorBeta,
    NoPriorSigma,
    compute_ss_from_pre,
    KeepReflection,
    DoKeepReflection,
    GeneralSigma,
    standardize_reg,
    SizeAndShapeMCMC_p2withreflection,
    Valuep2,
    compute_designmatrix,
    compute_dessignmatrix


end # module
