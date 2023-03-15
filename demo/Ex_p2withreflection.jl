### ### ### ### ### 
### packages
### ### ### ### ### 

using Random, Distributions, LinearAlgebra, StatsBase
using DataFrames, StatsModels, CategoricalArrays
using BayesSizeAndShape

### ### ### ### ### 
### simulations
### ### ### ### ### 

n::Int64 = 100;
p::Int64 = 2;
k::Int64 = 10;
d::Int64 = 3;

# regression
reg::Matrix{Float64} = zeros(Float64, k*d, p);
reg[:] = rand(Normal(0.0, 5.0), prod(size(reg)));
standardize_reg(reg::Matrix{Float64}, Valuep2())


zmat = DataFrame(
    x1 = rand(Normal(10.0,1.0 ),n),
    x2 = sample(["A", "B"],n)
)
zmat[:,1] = (zmat[:,1] .- mean(zmat[:,1])) ./ std(zmat[:,1])
zmat.x2 = categorical(zmat.x2)
zmat_modmat_ModelFrame = ModelFrame(@formula(1 ~ 1+x1+x2), zmat);
zmat_modmat = ModelMatrix(zmat_modmat_ModelFrame).m
design_matrix = compute_designmatrix(zmat_modmat, k); # dimensions  k, k * d, n


# covariance
sigma::Symmetric{Float64,Matrix{Float64}} = Symmetric(rand(InverseWishart(k + 2, 5.0 * Matrix{Float64}(I, k, k))));


dataset_complete = zeros(Float64,k,p,n);
dataset = zeros(Float64, k, p, n);
for i_n = 1:n
    for i_p = 1:p
        dataset_complete[:, i_p, i_n] = rand(MvNormal(design_matrix[:, :, i_n] * reg[:, i_p], sigma))
        
    end
end
rmat = compute_ss_from_pre(dataset_complete, dataset, true);

### ### ### ### ### 
### MCMC
### ### ### ### ### 
betaout, sigmaout = SizeAndShapeMCMC(;
    dataset = dataset,
    fm = @formula(1 ~ 1+ x1 + x2),
    covariates = zmat,
    iterations=(iter=1000, burnin=200, thin=2),
    betaprior = Normal(0.0, 10000.0),
    sigmaprior=InverseWishart(k + 2, 5.0 * Matrix{Float64}(I, k, k)),
    beta_init= zeros(Float64, k * d, p),
    sigma_init = Symmetric(Matrix{Float64}(I, k, k)),
    rmat_init = reshape(vcat([Matrix{Float64}(I, p, p)[:] for i = 1:n]...), (p, p, n)),
    reflection = true
);

