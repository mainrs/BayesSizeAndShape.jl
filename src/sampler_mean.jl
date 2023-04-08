#### #### #### #### #### #### 
#### BETA
#### #### #### #### #### #### 

function sampler_mean(

    datamodel::SSDataType{R,RL,VP,RS,IC},
    data_mcmc::MCMCData,
    mean_model::TypeModelMean,
    mean_mcmc::MCMCTypeModelMean,
    covariance_model::TypeModelCoVariance,
    covariance_mcmc::MCMCTypeModelCoVariance
    ) where {R<:Reflection,RL<:RemoveLocation,VP<:ValueP,RS<:RemoveSize,IC<:IdentifiabilityConstraint}

   
    error("sampler_mean")

end


function sampler_mean(
    datamodel::SSDataType{R,RL,VP,RS,IC},
    data_mcmc::MCMCData,
    mean_model::TypeModelMean,
    mean_mcmc::MCMCLinearMean{<:NoPriorBeta},
    covariance_model::TypeModelCoVariance,
    covariance_mcmc::MCMCTypeModelCoVariance
    ) where {R<:Reflection,RL<:RemoveLocation,VP<:ValueP,RS<:RemoveSize,IC<:IdentifiabilityConstraint}

    # No sample beta

end

function sampler_mean(
    datamodel::SSDataType{R,RL,VP,DoNotRemoveSize,IC},
    data_mcmc::MCMCNormalDataKeepSize{R,RL,VP},
    mean_model::LinearMean,
    mean_mcmc::MCMCLinearMean{<:Normal},
    covariance_model::GeneralCoVarianceIndependentDimension,
    covariance_mcmc::MCMCGeneralCoVarianceIndependentDimension
    ) where {R<:Reflection,VP<:ValueP,IC<:IdentifiabilityConstraint,RL<:DoRemoveLocation}
    

    betaMCMC = mean_mcmc.beta_mcmc
    invMat = covariance_mcmc.invcovariance_mcmc

    prior = mean_mcmc.prior_beta
    designmatrix = mean_model.designmatrix
    yrdata = data_mcmc.yr_mcmc

    kd::Int64 = size(betaMCMC,1)
    p::Int64  = datamodel.p
    n::Int64  = datamodel.n

    

    Vp::Matrix{Float64} = zeros(Float64, kd, kd)
    Mp::Vector{Float64} = zeros(Float64, kd)

    for ip = 1:p

        Vp[:, :] = Diagonal([1.0 / params(prior)[2]^2 for i = 1:kd])
        Mp[:] .= params(prior)[1] / params(prior)[2]^2

        for j = 1:n
            Vp[:, :] += transpose(designmatrix[:, :, j]) * invMat * designmatrix[:, :, j]
            Mp[:, :] += transpose(designmatrix[:, :, j]) * invMat * yrdata[:,ip,j]
        end

        Vp = Symmetric(inv(Vp))
        Mp = Vp*Mp

        betaMCMC[:,ip] = rand(MvNormal(Mp,Vp))
       
    end
    for ip = 1:p

        for j = 1:n
            mean_mcmc.mean_mcmc[:, ip, j] = designmatrix[:, :, j] * betaMCMC[:, ip]
        end
    end

end


#function sampler_beta(xdata::Array{Float64,3}, designmatrix::Array{Float64,3}, prior::NoPriorBeta, betaMCMC::Matrix{Float64}, sigmaMCMC::Matrix{Float64}, valp::Valuep, betastandMCMC::Matrix{Float64})

#end
#function sampler_beta(xdata::Array{Float64,3}, designmatrix::Array{Float64,3}, prior::Normal, betaMCMC::Matrix{Float64}, sigmaMCMC::Matrix{Float64}, valp::Valuep, betastandMCMC::Matrix{Float64})

    
    
#    #standardize_reg(betastandMCMC, valp)

#end
