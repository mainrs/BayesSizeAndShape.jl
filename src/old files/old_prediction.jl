
function predictmean(;beta::Array{Float64,3}, 
    dataset::Array{Float64,3}, 
    fm::FormulaTerm = @formula(1~ 1),
    covariates::DataFrame)
    
    k::Int64 = size(dataset, 1)
    n::Int64 = size(dataset, 3)
    p::Int64 = size(dataset, 2)

    mean_pred = [zeros(Float64, size(beta,1), k,p) for i = 1:n];

    
    covariates_copy = deepcopy(covariates)
    for i = 1:size(covariates_copy,2)
        if  isa(covariates_copy[:,i], CategoricalArray)

        elseif isa(covariates_copy[1,i], Real)
            covariates_copy[:,i] = (covariates_copy[:,i] .- mean(covariates_copy[:,i])) ./ std(covariates_copy[:,i])
        else
        
            error("Only factors or Real variables are allowed in covariates")
        end
    end

    designmatrix_v2_app = ModelFrame(fm, covariates_copy);
    designmatrix_v2 = ModelMatrix(designmatrix_v2_app).m
    #colnames_modelmatrix = coefnames(designmatrix_v2_app)

    #for i = 2:size(designmatrix_v2, 2)

    #    designmatrix_v2[:, i] = designmatrix_v2[:, i] .- mean(designmatrix_v2[:, i])

    #end
    designmatrix = compute_designmatrix(designmatrix_v2, k) # dimensions  k, k * d, n

    
    for select_obs = 1:n

        for i = 1:size(beta,1)

            mean_pred[select_obs][i,:,:] = designmatrix[:,:,select_obs]*beta[i,:,:]
    
        end
    end
    
    return mean_pred

end



function predictmean(modeloutput::SizeAndShapeModelOutput{KeepReflection, GeneralSigma, Valuep2, DB, DS, FF}) where {DB<:ContinuousUnivariateDistribution,DS<:ContinuousMatrixDistribution, FF<:FormulaTerm}
    

    covariates = modeloutput.covariates.covariates
    k::Int64 = modeloutput.indices.k
    n::Int64 = size(covariates,1)
    p::Int64 = modeloutput.indices.p

    beta = modeloutput.mcmcoutputArrays.beta
    rmat = modeloutput.mcmcoutputArrays.rmat

    mean_pred = [zeros(Float64, size(beta,1), k,p) for i = 1:n];

    
    #covariates_copy = deepcopy(covariates)
    #for i = 1:size(covariates_copy,2)
    #    if  isa(covariates_copy[:,i], CategoricalArray)

    #    elseif isa(covariates_copy[1,i], Real)
    #        #covariates_copy[:,i] = (covariates_copy[:,i] .- mean(covariates_copy[:,i])) ./ std(covariates_copy[:,i])
    #    else
        
    #        error("Only factors or Real variables are allowed in covariates")
    #    end
    #end

    #designmatrix_v2_app = ModelFrame(modeloutput.covariates.fm, covariates_copy);
    #designmatrix_v2 = ModelMatrix(designmatrix_v2_app).m
    


    #designmatrix = compute_designmatrix(designmatrix_v2, k) # dimensions  k, k * d, n

    designmatrix, d, colnames_modelmatrix, designmatrix_v2  = create_designmatrix(modeloutput.covariates.fm, covariates, k)
    for select_obs = 1:n

        for i = 1:size(beta,1)

            mean_pred[select_obs][i,:,:] = designmatrix[:,:,select_obs]*beta[i,:,:]*rmat[i,:,:,select_obs]
    
        end
    end
    
    return mean_pred

end