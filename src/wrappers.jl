
"""
    SizeAndShapeMCMC_p2withreflection(;
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
        rmat_init::Array{Float64,3}
)
Posterior sample from the size-and-shape model for two-dimensional data with reflection information. 


# Arguments
Let 
* n be the number of shapes;
* k be the number of landmarks (for the pre-form matrix **X**);
* p be the dimension of each landmark (only p=2 is implemented)
* d be the number of covariates to use;


The arguments are

- `dataset::Array{Float64,3}`: Array with the data - dimension (k,p,n) of size-and-shape data. Use the function `compute_ss_from_pre` to obtain the size-and-shape data from pre-forms
- `fm::FormulaTerm = @formula(1~ 1)`: a formula that specifies the model - the left-and size should be 1
- `covariates::DataFrame`: a DataFrame containing the covariates - dimension (n,d). The names used in `fm` must be column names of  `covariates`
- `iterations::NamedTuple{(:iter, :burnin, :thin),Tuple{Int64,Int64,Int64}}`: values of the iterations (iter), thin and burnin of the MCMC algorithm
- `betaprior::ContinuousUnivariateDistribution`: The prior on the regressive coefficients - only a Normal distribution is allowed
- `sigmaprior::ContinuousMatrixDistribution`: The prior on the covariance matrix - only an Inverse Wishart is allowed
- `beta_init::Matrix{Float64}`: initial values for the regressive coefficients - dimension (k*d,p) 
- `sigma_init::Symmetric{Float64,Matrix{Float64}}`: initial values for the covariance matrix - dimension (k,k) (it must be a valid covariance matrix) 
- `rmat_init::Array{Float64,3}`: initial values for the rotation matrices - dimension (p,p,n) (each [1:p, 1:p, i], i = 1,2,...,n, must be a valid rotation matrix) 
"""
function SizeAndShapeMCMC_p2withreflection(;
    dataset::Array{Float64,3}, 
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
    rmat_init::Array{Float64,3}
)

return mcmc(; 
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

end