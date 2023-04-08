#### #### #### #### #### #### 
#### BETA
#### #### #### #### #### #### 

function sampler_latentobservations(
    iterMCMC::Int64,
    datamodel::SSDataType{R,RL,VP,RS,IC},
    data_mcmc::MCMCData,
    mean_model::TypeModelMean,
    mean_mcmc::MCMCTypeModelMean,
    covariance_model::TypeModelCoVariance,
    covariance_mcmc::MCMCTypeModelCoVariance
    ) where {R<:Reflection,RL<:RemoveLocation,VP<:ValueP,RS<:RemoveSize,IC<:IdentifiabilityConstraint}

   
    error("sampler_latentobservations")

end


function sampler_latentobservations(
    iterMCMC::Int64,
    datamodel::SSDataType{R,RL,VP,RS,IC},
    data_mcmc::MCMCNormalDataKeepSize{R,RL,VP,<:DoNotSampleRmat},
    mean_model::TypeModelMean,
    mean_mcmc::MCMCTypeModelMean,
    covariance_model::TypeModelCoVariance,
    covariance_mcmc::MCMCTypeModelCoVariance
    ) where {R<:Reflection,RL<:RemoveLocation,VP<:ValueP,RS<:RemoveSize,IC<:IdentifiabilityConstraint}

    # No sample beta

end

function sampler_latentobservations(
    iterMCMC::Int64,
    datamodel::SSDataType{KeepReflection,RL,ValueP2,RS,IC},
    data_mcmc::MCMCNormalDataKeepSize{KeepReflection,RL,ValueP2,<:DoSampleRmat},
    mean_model::TypeModelMean,
    mean_mcmc::MCMCTypeModelMean,
    covariance_model::TypeModelCoVariance,
    covariance_mcmc::MCMCTypeModelCoVariance
    ) where {RL<:RemoveLocation,RS<:RemoveSize,IC<:IdentifiabilityConstraint}



    #meanMCMC = mean_mcmc.mean_mcmc
    #sigmaMCMC = covariance_mcmc.covariance_mcmc
    #invMat = covariance_mcmc.invcovariance_mcmc
    
    

    #prior = covariance_mcmc.prior_covariance
    
    #yrdata = data_mcmc.yr_mcmc

    #betaMCMC = mean_mcmc.beta_mcmc
    invMat = covariance_mcmc.invcovariance_mcmc
    

    p::Int64  = datamodel.p
    n::Int64  = datamodel.n
    k::Int64  = datamodel.k


    for j = 1:n
    
        #mui = designmatrix[:, :, j] * betaMCMC[:, :]
    
        #@toggled_assert size(mui) == (k,p)
        #println(size(datamodel.ssdata[:,:,j]))
        Ai = transpose(mean_mcmc.mean_mcmc[:,:,j]) * invMat * datamodel.ssdata[:,:,j]
        #println(size(Ai))

        x1 = Ai[1, 1] + Ai[2, 2]
        x2 = Ai[2, 1] - Ai[1, 2]
        #x2 = Ai[1, 2] - Ai[2, 1]

        kvonmises = sqrt(x1^2 + x2^2)

        muvonmises = atan(x2 / kvonmises, x1 / kvonmises)
        #muvonmises = atan(x1 / kvonmises, x2 / kvonmises)

        data_mcmc.angles_mcmc[1, j] = rand(VonMises(muvonmises, kvonmises))

        data_mcmc.rmat_mcmc[1, 1, j] = cos(data_mcmc.angles_mcmc[1, j])
        data_mcmc.rmat_mcmc[1, 2, j] = -sin(data_mcmc.angles_mcmc[1, j])
        data_mcmc.rmat_mcmc[2, 1, j] = sin(data_mcmc.angles_mcmc[1, j])
        data_mcmc.rmat_mcmc[2, 2, j] = cos(data_mcmc.angles_mcmc[1, j])


        #xdata[:, :, j] = ydata[:, :, j] * transpose(rmatMCMC[:,:,j])
    end
    compute_yrdata(data_mcmc.yr_mcmc, datamodel.ssdata, data_mcmc.rmat_mcmc)

end



function sampler_latentobservations(
    iterMCMC::Int64,
    datamodel::SSDataType{R,RL,VP,RS,IC},
    data_mcmc::MCMCNormalDataKeepSize{R,RL,VP,<:DoNotSampleRmat},
    mean_model::TypeModelMean,
    mean_mcmc::MCMCTypeModelMean,
    covariance_model::TypeModelCoVariance,
    covariance_mcmc::MCMCTypeModelCoVariance
    ) where {R<:Reflection,RL<:RemoveLocation,VP<:ValueP,RS<:RemoveSize,IC<:IdentifiabilityConstraint}

    # No sample beta

end

function sampler_latentobservations(
    iterMCMC::Int64,
    datamodel::SSDataType{KeepReflection,RL,ValueP3,RS,IC},
    data_mcmc::MCMCNormalDataKeepSize{KeepReflection,RL,ValueP3,<:DoSampleRmat},
    mean_model::TypeModelMean,
    mean_mcmc::MCMCTypeModelMean,
    covariance_model::TypeModelCoVariance,
    covariance_mcmc::MCMCTypeModelCoVariance
    ) where {RL<:RemoveLocation,RS<:RemoveSize,IC<:IdentifiabilityConstraint}

    n::Int64  = datamodel.n
    for j = 1:n

        sampler_latentobservations(j, iterMCMC, datamodel, data_mcmc, mean_model,  mean_mcmc, covariance_model, covariance_mcmc)

    end
    compute_yrdata(data_mcmc.yr_mcmc, datamodel.ssdata, data_mcmc.rmat_mcmc)


end
function sampler_latentobservations(j::Int64,
    iterMCMC::Int64,
    datamodel::SSDataType{KeepReflection,RL,ValueP3,RS,IC},
    data_mcmc::MCMCNormalDataKeepSize{KeepReflection,RL,ValueP3,<:DoSampleRmat},
    mean_model::TypeModelMean,
    mean_mcmc::MCMCTypeModelMean,
    covariance_model::TypeModelCoVariance,
    covariance_mcmc::MCMCTypeModelCoVariance
    ) where {RL<:RemoveLocation,RS<:RemoveSize,IC<:IdentifiabilityConstraint}


    #### https://arxiv.org/abs/math/0503712
    #### Bayesian alignment using hierarchical models, with applications in protein bioinformatics


    

    p::Int64  = datamodel.p
    n::Int64  = datamodel.n
    k::Int64  = datamodel.k
    angles = data_mcmc.angles_mcmc
    invMat = covariance_mcmc.invcovariance_mcmc

    MHratio::Float64 = 0.0
    molt::Float64 = 0.0
    alpha_mean::Float64 = 0.0

    Ai = transpose(mean_mcmc.mean_mcmc[:,:,j]) * invMat * datamodel.ssdata[:,:,j]
        
    #### angle1 (-pi, pi) 
    x1 = (Ai[2,2] - sin(angles[2,j]) * Ai[1,3]) * cos(angles[3,j])
    x1 += (-Ai[2,3]- sin(angles[2,j])*Ai[1,2]) * sin(angles[3,j])
    x1 += cos(angles[2,j])*Ai[1,1]
    
    x2 = (-sin(angles[2,j])*Ai[2,3]-Ai[1,2]) * cos(angles[3,j])
    x2 += (Ai[1,3]- sin(angles[2,j])*Ai[2,2]) * sin(angles[3,j])
    x2 += cos(angles[2,j])*Ai[2,1]

    kvonmises = sqrt(x1^2 + x2^2)
    muvonmises = atan(x2 / kvonmises, x1 / kvonmises)
    angles[1, j] = rand(VonMises(muvonmises, kvonmises))

    #### angle2 (-pi/2, pi/2)
    x1 = sin(angles[1, j])*Ai[2,1] + cos(angles[1, j])*Ai[1,1]
    x1 += sin(angles[3, j])*Ai[3,2] + cos(angles[3, j])*Ai[3,3]

    x2 = (-sin(angles[3, j])*Ai[1,2]-cos(angles[3, j])*Ai[1,3])*cos(angles[1, j])
    x2 += (-sin(angles[3, j])*Ai[2,2]-cos(angles[3, j])*Ai[2,3])*sin(angles[1, j])
    x2 += Ai[3,1]

    

    # prop = rand(Uniform(angles[2, j]- lim, angles[2, j] + lim))
    


    prop = rand(Normal(angles[2, j],data_mcmc.sdprop_adapt[2,j]))

    if (prop> - (pi/Float64(2.0))) & (prop < (pi/Float64(2.0)))

        MHratio = x1*cos(prop) + x2*sin(prop) + log(cos(prop))
        MHratio -= x1*cos(angles[2, j]) + x2*sin(angles[2, j]) + log(cos(angles[2, j]))

        if rand(Uniform(0.0,1.0)) < exp(MHratio)
            angles[2, j] = prop
            #println("Acc")
        else
            #println("Rej")
        end
        data_mcmc.sumalpha[2,j] += min(exp(MHratio), 1.0)
    else
        data_mcmc.sumalpha[2,j] += 0.0
    end

    

    

    if (iterMCMC>data_mcmc.init_adapt[2,j]) & (iterMCMC<data_mcmc.end_adapt[2,j])
        
        molt = data_mcmc.a_adapt[2,j]/(data_mcmc.b_adapt[2,j] + iterMCMC)
        alpha_mean = data_mcmc.sumalpha[2,j]/data_mcmc.iter_adapt[2,j]

        data_mcmc.sdprop_adapt[2,j] = exp( log(data_mcmc.sdprop_adapt[2,j]) + molt*(alpha_mean - data_mcmc.accratio_adapt[2,j]) )

    end

    if mod(iterMCMC, data_mcmc.iter_adapt[2,j]) == 0

        data_mcmc.sumalpha[2,j] = 0.0

    end


    #### angle3 (-pi, pi) 
    x1 = (Ai[2,2] - sin(angles[2, j])*Ai[1,3])*cos(angles[1, j])
    x1 += (-sin(angles[2, j])*Ai[2,3]-Ai[1,2])*sin(angles[1, j])
    x1 += cos(angles[2,j])*Ai[3,3]

    x2 = (-Ai[2,3]-sin(angles[2, j])*Ai[1,2])*cos(angles[1, j])
    x2 += (Ai[1,3]-sin(angles[2, j])*Ai[2,2])*sin(angles[1, j])
    x2 += cos(angles[2,j])*Ai[3,2]

    kvonmises = sqrt(x1^2 + x2^2)
    muvonmises = atan(x2 / kvonmises, x1 / kvonmises)
    angles[3, j] = rand(VonMises(muvonmises, kvonmises))


    compute_rtridim_from_angles(j, data_mcmc.rmat_mcmc, angles)
    #data_mcmc.rmat_mcmc[1, 1, j] = cos(data_mcmc.angles_mcmc[1, j])
    #data_mcmc.rmat_mcmc[1, 2, j] = -sin(data_mcmc.angles_mcmc[1, j])
    #data_mcmc.rmat_mcmc[2, 1, j] = sin(data_mcmc.angles_mcmc[1, j])
    #data_mcmc.rmat_mcmc[2, 2, j] = cos(data_mcmc.angles_mcmc[1, j])


    #xdata[:, :, j] = ydata[:, :, j] * transpose(rmatMCMC[:,:,j])
    

end


function compute_rtridim_from_angles(i::Int64, rmat::Array{Float64,3},  angle::Matrix{Float64})



    M1::Matrix{Float64} = zeros(Float64,3,3)
    M2::Matrix{Float64} = zeros(Float64,3,3)
    M3::Matrix{Float64} = zeros(Float64,3,3)

    M1[1,1] = cos(angle[1,i])
    M1[2,2] = cos(angle[1,i])
    M1[1,2] = -sin(angle[1,i])
    M1[2,1] = sin(angle[1,i])
    M1[3,3] = Float64(1.0)

    M2[1,1] = cos(angle[2,i])
    M2[3,3] = cos(angle[2,i])
    M2[1,3] = -sin(angle[2,i])
    M2[3,1] = sin(angle[2,i])
    M2[2,2] = Float64(1.0)

    M3[2,2] = cos(angle[3,i])
    M3[3,3] = cos(angle[3,i])
    M3[2,3] = -sin(angle[3,i])
    M3[3,2] = sin(angle[3,i])
    M3[1,1] = Float64(1.0)
    
    rmat[:,:,i] = M1*M2*M3

    return nothing
end

#function sampler_beta(xdata::Array{Float64,3}, designmatrix::Array{Float64,3}, prior::NoPriorBeta, betaMCMC::Matrix{Float64}, sigmaMCMC::Matrix{Float64}, valp::Valuep, betastandMCMC::Matrix{Float64})

#end
#function sampler_beta(xdata::Array{Float64,3}, designmatrix::Array{Float64,3}, prior::Normal, betaMCMC::Matrix{Float64}, sigmaMCMC::Matrix{Float64}, valp::Valuep, betastandMCMC::Matrix{Float64})

    
    
#    #standardize_reg(betastandMCMC, valp)

#end
