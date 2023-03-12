
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
        rmat_init::Array{Float64,3}
)
Posterior sample form the size-and-shape model for a two-dimensional data with reflection information. 


# Arguments
Let 
* n be the number of shapes;
* k be the number of landmark (for the pre-form matrix **X**);
* d be the dimension of each landmark (only p=2 is implemented)
* p be the number of covariates to use;


The arguments are

- `dataset::Array{Float64,3}`: Array with the data - dimension (k,p,n) 
- `fm::FormulaTerm = @formula(1~ 1)`: a formula that specifies the model -  the left-and size  shoudl be 1
- `covariates::DataFrame`: a DataFrame containing the covariates - dimension (n,d)
- `iterations::NamedTuple{(:iter, :burnin, :thin),Tuple{Int64,Int64,Int64}}`: values of the iterations (iter), thin and burnin of the MCMC algorithm
- `betaprior::ContinuousUnivariateDistribution`: The prior on the regressive coefficients -  only a Normal distribution is allowed
- `sigmaprior::ContinuousMatrixDistribution`: The prior on the covariance matrix -  only an Inverse Wishart is allowed
- `beta_init::Matrix{Float64}`: initial values for the regressive coefficients -  dimension (k*d,p) 
- `sigma_init::Symmetric{Float64,Matrix{Float64}}`: initial values for the covariance matrix -  dimension (k,k) (it must be a valid covariance matrix) 
- `rmat_init::Array{Float64,3}`: initial values for the rotation matrices -  dimension (p,p,n) (each [:,:,i] matrix must be a valid rotation matrix) 
...
"""
function SizeAndShapeMCMC(;
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