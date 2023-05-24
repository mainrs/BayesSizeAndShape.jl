function compute_dessignmatrix(zmat::Matrix{Float64}, k::Int64)::Array{Float64,3}
    println("compute_dessignmatrix is deprecated. Please use compute_designmatrix instead")
    return compute_designmatrix(zmat, k)

end

function compute_dessignmatrix(zmat::DataFrame, k::Int64)::Array{Float64,3}
    println("compute_dessignmatrix is deprecated. Please use compute_designmatrix instead")
    return compute_designmatrix(zmat, k)

end

function compute_ss_from_pre(xdata::Array{Float64,3},ydata::Array{Float64,3}, keep_reflection::Bool)::Array{Float64,3}

    println("compute_ss_from_pre is deprecated. Please use compute_ss_helmertz! instead")

    return compute_ss!(xdata,ydata, keep_reflection)

end


"""
    SizeAndShapeMCMC(;
        dataset::Array{Float64,3}, 
        fm::FormulaTerm = @formula(1~ 1),
        covariates::DataFrame,
        iterations::NamedTuple{(:iter, :burnin, :thin),Tuple{Int64,Int64,Int64}} =(
            iter=1000,
            burnin=200,
            thin=2
        ),
        betaprior::ContinuousUnivariateDistribution = Normal(0.0, 10000.0),
        sigmaprior::ContinuousMatrixDistribution,
        beta_init::Matrix{Float64},
        sigma_init::Symmetric{Float64,Matrix{Float64}},
        rmat_init::Array{Float64,3},
        reflection::Bool = true
        )

Posterior samples from the size-and-shape model - in this version, only two-dimensional data with reflection information are allowed. 
The functions returns 2 Dataframes, with the posterior samples of the regressive coefficients and elements of the covariance matrix.


# Arguments
Let 
* n be the number of shapes;
* k be the number of landmarks (on the pre-form matrix);
* p be the dimension of each landmark (only p=2 is implemented)
* d be the number of covariates to use + 1 (or number of regressive coefficients for each landmark-dimension, intercept included);


The arguments are

- `dataset::Array{Float64,3}`: Array with the data - dimension (k,p,n) of size-and-shape data. Use the function `compute_ss_from_pre` to obtain the size-and-shape data from pre-forms
- `fm::FormulaTerm = @formula(1~ 1)`: a formula that specifies the model - the left-and size should be 1
- `covariates::DataFrame`: a DataFrame containing the covariates - dimension (n,d). The names used in `fm` must be column names of  `covariates`. Only Real And Factor covariates are allowed. The numeric column are standardized internally.
- `iterations::NamedTuple{(:iter, :burnin, :thin),Tuple{Int64,Int64,Int64}}`: values of the iterations (iter), thin and burnin of the MCMC algorithm
- `betaprior::ContinuousUnivariateDistribution`: The prior on the regressive coefficients - only a Normal distribution is allowed
- `sigmaprior::ContinuousMatrixDistribution`: The prior on the covariance matrix - only an Inverse Wishart is allowed
- `beta_init::Matrix{Float64}`: initial values for the regressive coefficients - dimension (k*d,p) 
- `sigma_init::Symmetric{Float64,Matrix{Float64}}`: initial values for the covariance matrix - dimension (k,k) (it must be a valid covariance matrix) 
- `rmat_init::Array{Float64,3}`: initial values for the rotation matrices - dimension (p,p,n) (each [1:p, 1:p, i], i = 1,2,...,n, must be a valid rotation matrix) 
- `reflection::Bool`: true for a model with reflection information.
"""
function SizeAndShapeMCMC(;
    dataset_init::Array{Float64,3}, 
    fm::FormulaTerm = @formula(1~ 1),
    covariates::DataFrame, #
    iterations::NamedTuple{(:iter, :burnin, :thin),Tuple{Int64,Int64,Int64}} =(
        iter=1000,
        burnin=200,
        thin=2
    ),
    betaprior::ContinuousUnivariateDistribution = Normal(0.0, 10000.0),
    sigmaprior::ContinuousMatrixDistribution,
    beta_init::Matrix{Float64},
    sigma_init::Symmetric{Float64,Matrix{Float64}},
    rmat_init::Array{Float64,3},
    reflection::Bool = true,
    is_data_sizeandshape::Bool = true
)

    println("SizeAndShapeMCMC is deprecated. Please use generalSizeAndShapeMCMC")
    p::Int64 = size(dataset, 2)
    if p != 2
        error("the current implementation allows only `size(dataset, 2) = 2 (2-dimensional data)`")
    end
    if reflection == false
        error("the current implementation allows only `reflection = true`")
    else
         mcmcout =generalSizeAndShapeMCMC(; 
            dataset = dataset,
            fm = fm,
            covariates = covariates,
            iterations = iterations,
            betaprior = betaprior,
            sigmaprior = sigmaprior,
            beta_init = beta_init,
            sigma_init = sigma_init,
            rmat_init = rmat_init,
            dormat = true,
            reflection = KeepReflection(),
            sigmatype = GeneralSigma()
        )
        
        return mcmcout.mcmcoutputArrays.beta, mcmcout.mcmcoutputArrays.sigma
    end
    

end