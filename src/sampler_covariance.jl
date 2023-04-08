#### #### #### #### #### #### 
#### BETA
#### #### #### #### #### #### 

function sampler_covariance(

    datamodel::SSDataType{R,RL,VP,RS,IC},
    data_mcmc::MCMCData,
    mean_model::TypeModelMean,
    mean_mcmc::MCMCTypeModelMean,
    covariance_model::TypeModelCoVariance,
    covariance_mcmc::MCMCTypeModelCoVariance
    ) where {R<:Reflection,RL<:RemoveLocation,VP<:ValueP,RS<:RemoveSize,IC<:IdentifiabilityConstraint}

   
    error("sampler_covariance")

end


function sampler_covariance(
    datamodel::SSDataType{R,RL,VP,RS,IC},
    data_mcmc::MCMCData,
    mean_model::TypeModelMean,
    mean_mcmc::MCMCLinearMean,
    covariance_model::TypeModelCoVariance,
    covariance_mcmc::MCMCGeneralCoVarianceIndependentDimension{<:NoPriorSigma}
    ) where {R<:Reflection,RL<:RemoveLocation,VP<:ValueP,RS<:RemoveSize,IC<:IdentifiabilityConstraint}

    # No sample beta
    #println("No Prir")
end

function sampler_covariance(
    datamodel::SSDataType{R,RL,VP,DoNotRemoveSize,IC},
    data_mcmc::MCMCNormalDataKeepSize{R,RL,VP},
    mean_model::LinearMean,
    mean_mcmc::MCMCLinearMean,
    covariance_model::GeneralCoVarianceIndependentDimension,
    covariance_mcmc::MCMCGeneralCoVarianceIndependentDimension{<:InverseWishart}
    ) where {R<:Reflection,VP<:ValueP,IC<:IdentifiabilityConstraint,RL<:DoRemoveLocation}

    

    meanMCMC = mean_mcmc.mean_mcmc
    sigmaMCMC = covariance_mcmc.covariance_mcmc
    invMat = covariance_mcmc.invcovariance_mcmc
    
    

    prior = covariance_mcmc.prior_covariance
    
    yrdata = data_mcmc.yr_mcmc

    #kd::Int64 = size(betaMCMC,1)
    p::Int64  = datamodel.p
    n::Int64  = datamodel.n




    nup::Float64 = params(prior)[1] + Float64(n*p) 
    Psip = deepcopy(params(prior)[2].mat)
    app = yrdata[:, :,:] - meanMCMC[:,:,:]
    for ip = 1:p
        for j = 1:n
            
            Psip[:,:] += app[:,ip,j]*transpose(app[:,ip,j])
            
        end
    end
    
    #println("a1= ", yrdata[:,:,1])
    #println("a2 = ", meanMCMC[:,:,1])
    #println("ssdata = ", data_mcmc.ssdata[:,:,1])
    

    
    #println("")
    #println([n,p])
    #println("")
    #println(nup)
    #println(Psip)



    sigmaMCMC.data[:, :] = rand(InverseWishart(nup, Symmetric(Psip).data))

#println(sigmaMCMC.data)

    cc = cholesky(sigmaMCMC)
    covariance_mcmc.invcovariance_mcmc.data[:,:] = inv(cc)[:,:]
    covariance_mcmc.logdeterminant_mcmc[:]= [log(det(cc))]

    #error("")
    return nothing
    
end


#function sampler_beta(xdata::Array{Float64,3}, designmatrix::Array{Float64,3}, prior::NoPriorBeta, betaMCMC::Matrix{Float64}, sigmaMCMC::Matrix{Float64}, valp::Valuep, betastandMCMC::Matrix{Float64})

#end
#function sampler_beta(xdata::Array{Float64,3}, designmatrix::Array{Float64,3}, prior::Normal, betaMCMC::Matrix{Float64}, sigmaMCMC::Matrix{Float64}, valp::Valuep, betastandMCMC::Matrix{Float64})

    
    
#    #standardize_reg(betastandMCMC, valp)

#end
