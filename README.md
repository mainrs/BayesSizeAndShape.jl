# BayesSizeAndShape.jl

<!--![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)<!--
![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-stable-green.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-retired-orange.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-archived-red.svg)
![Lifecycle](https://img.shields.io/badge/lifecycle-dormant-blue.svg) -->
<!--[![Build Status](https://travis-ci.com/GianlucaMastrantonio/BayesSizeAndShape.jl.svg?branch=master)](https://travis-ci.com/GianlucaMastrantonio/BayesSizeAndShape.jl)
[![codecov.io](http://codecov.io/github/GianlucaMastrantonio/BayesSizeAndShape.jl/coverage.svg?branch=master)](http://codecov.io/github/GianlucaMastrantonio/BayesSizeAndShape.jl?branch=master)-->
<!--
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://GianlucaMastrantonio.github.io/BayesSizeAndShape.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-master-blue.svg)](https://GianlucaMastrantonio.github.io/BayesSizeAndShape.jl/dev)
-->

To install the package, simply do
```julia
julia> ]

pkg> add BayesSizeAndShape
```
at the julia prompt.

# **BayesSizeAndShape**

This package implements a Bayesian regression for size-and-shape data, based on the CITE.
At the present moment, the package implements a model for two-dimensional data with reflection information. The function to use is `SizeAndShapeMCMC`

### **Basic Usage**

In the **demo** directory, there is a julia file with an example of how to implement the model, which we will describe also here. 

Let 
* n be the number of shapes;
* k be the number of landmark (for the pre-form matrix **X**);
* d be the dimension of each landmark (only p=2 is implemented)
* p be the number of covariates to use;

We first simulate the regressive coefficients, which are in the `reg` object -  `reg` must be of  dimension (kd, p):
```julia
n::Int64 = 100;
p::Int64 = 2;
k::Int64 = 10;
d::Int64 = 3;

# regression
reg::Matrix{Float64} = zeros(Float64, k*d, p);
reg[:] = rand(Normal(0.0, 5.0), prod(size(reg)));
```
The regressive coefficients are standardized with the function `standardize_reg`. The second argument is here used to specify the dimension of each landmark, `Value2()` is used for two-dimensional data, while `Valuep3()` (not yet implemented) for three-dimensional data:  
```julia
standardize_reg(reg::Matrix{Float64}, Valuep2())
```

A design matrix is simulated and the function `compute_dessignmatrix` is used to transform to be used in the model - `zmat`  must be of  dimension (n, d):
```julia
zmat = zeros(Float64,  n,d);
zmat[:] = rand(Normal(0.0, 5.0), prod(size(zmat)));
zmat[:,1] .= 1.0;
zmat = DataFrame(zmat, :auto)
for i = 2:size(zmat,2)

    zmat[:, i] = zmat[:, i] .- mean(zmat[:, i])

end
design_matrix = compute_dessignmatrix(zmat, k);
```

The pre-form matrix **X**, here saved in the object `dataset_complete`, is simulated from a multivariate normal, and its size-and-shape version, contained in the object `dataset`, is obtained by using the function `compute_ss_from_pre`. The third argument of the function is used to specify if reflection must be kept (`true` is the only available option in this implementation)
```julia

sigma::Symmetric{Float64,Matrix{Float64}} = Symmetric(rand(InverseWishart(k + 2, 5.0 * Matrix{Float64}(I, k, k))));

dataset_complete = zeros(Float64,k,p,n);
dataset = zeros(Float64, k, p, n);
for i_n = 1:n
    for i_p = 1:p
        dataset_complete[:, i_p, i_n] = rand(MvNormal(design_matrix[:, :, i_n] * reg[:, i_p], sigma))
        
    end
end
rmat = compute_ss_from_pre(dataset_complete, dataset, true);
```

Posterior samples are obtained with the function `SizeAndShapeMCMC`. To specify the regressive formula, you can use `@formula`, where on the left-hand size there must be 1 and in the right-hand size is the actual regressive formula.
```julia
betaOUT, sigmaOUT, rmatOUT, angleOUT = SizeAndShapeMCMC(;
    dataset = dataset,
    fm = @formula(1 ~ x1+x2+x3),
    covariates = zmat,
    iterations = (iter=1000, burnin=200, thin=2),
    betaprior = Normal(0.0, 10000.0),
    sigmaprior = InverseWishart(k + 2, 5.0 * Matrix{Float64}(I, k, k)),
    beta_init = zeros(Float64, k * d, p),
    sigma_init = Symmetric(Matrix{Float64}(I, k, k)),
    rmat_init = reshape(vcat([Matrix{Float64}(I, p, p)[:] for i = 1:n]...), (p, p, n))
);
```
The objects `betaOUT`, `sigmaOUT`, `rmatOUT`, `angleOUT` contains the posterior samples (the samples are on the first indices of these objects). `angleOUT` contains the angles used to compute the matrix **R** which are in `rmatOUT`.

## Citing

See `CITATION.bib`