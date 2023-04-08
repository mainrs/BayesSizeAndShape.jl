#### #### #### #### #### #### 
#### Model OUT
#### #### #### #### #### #### 
abstract type outputMCMCSizeAndShape end

struct SizeAndShapeModelOutput{R<:Reflection,RL<:RemoveLocation,VP<:ValueP,RS<:RemoveSize,IC<:IdentifiabilityConstraint
    ,DM<:MCMCData, MT<:TypeModelMean,MM<:MCMCTypeModelMean,CT<:TypeModelCoVariance,CM<:MCMCTypeModelCoVariance,PS<:MCMCObjectOUT
    } <: outputMCMCSizeAndShape

    

        datatype::SSDataType{R,RL,VP,RS,IC}
        datamcmc::DM
        meantype::MT
        meanmcmc::MM
        covtype::CT
        covmcmc::CM
        posteriorsamples::PS

        mcmciterations::NamedTuple{(:iter, :burnin, :thin, :savedsamples),Tuple{Int64,Int64,Int64, Int64}}

        function SizeAndShapeModelOutput(
            datatype::SSDataType{R,RL,VP,RS,IC}, 
            datamcmc::DM, 
            meantype::MT, 
            meanmcmc::MM, 
            covtype::CT, 
            covmcmc::CM, 
            posteriorsamples::PS,
            mcmciterations::NamedTuple{(:iter, :burnin, :thin, :savedsamples),Tuple{Int64,Int64,Int64, Int64}}
            ) where {R<:Reflection,RL<:RemoveLocation,VP<:ValueP,RS<:RemoveSize,IC<:IdentifiabilityConstraint
            ,DM<:MCMCData, MT<:TypeModelMean,MM<:MCMCTypeModelMean,CT<:TypeModelCoVariance,CM<:MCMCTypeModelCoVariance,PS<:MCMCObjectOUT
            }

            new{R,RL,VP,RS,IC,DM, MT,MM,CT,CM,PS}(datatype, datamcmc, meantype, meanmcmc, covtype, covmcmc,posteriorsamples, mcmciterations)
        end

end
#abstract type outputMCMCSizeAndShape end
#struct SizeAndShapeModelOutput{TR<:Reflection,TS<:SigmaType,TP<:Valuep,DB<:ContinuousUnivariateDistribution,DS<:ContinuousMatrixDistribution, FF<:FormulaTerm} <: outputMCMCSizeAndShape

#    mcmcoutput::NamedTuple{(:beta, :sigma, :rmat, :angle), Tuple{DataFrame, DataFrame, DataFrame, DataFrame}}

#    mcmcoutputArrays::NamedTuple{(:beta, :sigma, :rmat, :angle), Tuple{Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 4}, Array{Float64, 3}}}

#    covariates::NamedTuple{(:colnames_modelmatrix, :fm, :covariates, :designmatrix_step1, :designmatrix_step2), Tuple{Vector{String}, FF, DataFrame, Matrix{Float64}, Array{Float64, 3}}   }
    
#    dataset::Array{Float64,3}; 

#    modeltypes::NamedTuple{(:dormat, :reflection, :sigmatype, :betaprior, :sigmaprior, :valp), Tuple{Bool, TR, TS, DB,DS, TP}  }

#    iterations::NamedTuple{(:iter, :burnin, :thin, :savedsamples),Tuple{Int64,Int64,Int64, Int64}};
    
#    indices::NamedTuple{(:k, :p, :n, :d),Tuple{Int64,Int64,Int64, Int64}};
 
#    function SizeAndShapeModelOutput(
#        mcmcoutput::NamedTuple{(:beta, :sigma, :rmat, :angle), Tuple{DataFrame, DataFrame, DataFrame, DataFrame}},
#        mcmcoutputArrays::NamedTuple{(:beta, :sigma, :rmat, :angle), Tuple{Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 4}, Array{Float64, 3}}},
#        covariates::NamedTuple{(:colnames_modelmatrix, :fm, :covariates, :designmatrix_step1, :designmatrix_step2), Tuple{Vector{String}, FF, DataFrame, Matrix{Float64}, Array{Float64, 3}}   },
#        dataset::Array{Float64,3},
#        modeltypes::NamedTuple{(:dormat, :reflection, :sigmatype, :betaprior, :sigmaprior, :valp), Tuple{Bool, TR, TS, DB,DS, TP}  },
#        iterations::NamedTuple{(:iter, :burnin, :thin, :savedsamples),Tuple{Int64,Int64,Int64, Int64}},
#        indices::NamedTuple{(:k, :p, :n, :d),Tuple{Int64,Int64,Int64, Int64}}) where {TR<:Reflection,TS<:SigmaType,TP<:Valuep,DB<:ContinuousUnivariateDistribution,DS<:ContinuousMatrixDistribution,FF<:FormulaTerm} 

#        new{TR,TS,TP,DB,DS, FF}(mcmcoutput, mcmcoutputArrays, covariates, dataset, modeltypes, iterations, indices)
   
#    end
#end

#### #### #### #### #### #### 
#### FUNCTIONS
#### #### #### #### #### #### 

#function create_output(beta::Array{Float64, 3},
#    sigma::Array{Float64, 3},
#    rmat::Array{Float64, 4},
#    angle::Array{Float64, 3},
#    beta_nonid::Array{Float64, 3},
#    rmat_nonid::Array{Float64, 4},
#    angle_nonid::Array{Float64, 3},
#    colnames_modelmatrix::Vector{String},
#    fm::FormulaTerm,
#    betaprior::ContinuousUnivariateDistribution,
#    sigmaprior::ContinuousMatrixDistribution,
#    k::Int64,
#    p::Int64,
#    n::Int64,
#    d::Int64,
#    dormat::Bool,
#    reflection::Reflection,
#    sigmatype::SigmaType,
#    dataset::Array{Float64,3}, 
#    covariates::DataFrame,
#    iterations::NamedTuple{(:iter, :burnin, :thin),Tuple{Int64,Int64,Int64}},
#    designmatrix_step1::Array{Float64, 3},
#    designmatrix_step2::Matrix{Float64})

#    p2KeepReflectionGeneralSigma(beta::Array{Float64, 3},
#            sigma::Array{Float64, 3},
#            rmat::Array{Float64, 4},
#            angle::Array{Float64, 3},
#            beta_nonid::Array{Float64, 3},
#            rmat_nonid::Array{Float64, 4},
#            angle_nonid::Array{Float64, 3},
#            colnames_modelmatrix::Vector{String},
#            fm::FormulaTerm,
#            betaprior::ContinuousUnivariateDistribution,
#            sigmaprior::ContinuousMatrixDistribution,
#            k::Int64,
#            p::Int64,
#            n::Int64,
#            d::Int64,
#            dormat::Bool,
#            reflection::Reflection,
#            sigmatype::SigmaType,
#            dataset::Array{Float64,3}, 
#            covariates::DataFrame,
#            iterations::NamedTuple{(:iter, :burnin, :thin),Tuple{Int64,Int64,Int64}},
#            designmatrix_step1::Array{Float64, 3},
#            designmatrix_step2::Matrix{Float64})

#end