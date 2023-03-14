# **BayesSizeAndShape**

This package implements a Bayesian regression for size-and-shape data, based on "Di Noia, A., Mastrantonio, G., and Jona Lasinio, G., “Bayesian Size-and-Shape regression modelling”, <i>arXiv e-prints</i>, 2023".
At the present moment, the package implements a model for two-dimensional data with reflection information. The function to use is `SizeAndShapeMCMC_p2withreflection`

### **Basic Usage**

In the **demo** directory, there is a **julia** file **Ex_p2withreflection.jl** with an example of how to implement the model, which we will describe also here. 

Let 
* n be the number of shapes;
* k be the number of landmarks (for the pre-form matrix **X**);
* p be the dimension of each landmark (only p=2 is implemented)
* d be the number of covariates to use;

We first simulate the regressive coefficients, which are in the `reg` object - `reg` must be of dimension (kd, p):
```julia
n::Int64 = 100;
p::Int64 = 2;
k::Int64 = 10;
d::Int64 = 3;

# regression
reg::Matrix{Float64} = zeros(Float64, k*d, p);
reg[:] = rand(Normal(0.0, 5.0), prod(size(reg)));
```
The regressive coefficients are standardized with the function `standardize_reg`, for identifiability purposes. The second argument is here used to specify the dimension of each landmark, `Value2()` is used for two-dimensional data, while `Valuep3()` (not yet implemented) is for three-dimensional data:  
```julia
standardize_reg(reg::Matrix{Float64}, Valuep2())
```

A `DataFrame` named `zmat`, containing the covariates is simulated, and the function `compute_designmatrix` is used to compute the design matrix - `zmat` must be of dimension (n, d):
```julia
zmat = zeros(Float64,  n,d);
zmat[:] = rand(Normal(0.0, 5.0), prod(size(zmat)));
zmat[:,1] .= 1.0;
zmat = DataFrame(zmat, :auto)
for i = 2:size(zmat,2)

    zmat[:, i] = zmat[:, i] .- mean(zmat[:, i])

end
design_matrix = compute_designmatrix(zmat, k);
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

Posterior samples are obtained with the function `SizeAndShapeMCMC_p2withreflection`. To specify the regressive formula, you can use `@formula`, where on the left-hand side there must be 1 and on the right-hand side is the actual regressive formula.
```julia
betaOUT, sigmaOUT, rmatOUT, angleOUT = SizeAndShapeMCMC_p2withreflection(;
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
The objects `betaOUT`, `sigmaOUT`, `rmatOUT`, `angleOUT` contain the posterior samples. `angleOUT` contains the angles used to compute the matrix **R** in `rmatOUT`.

### **OUTPUT**

 - `betaOUT`: dimensions (l,kd, p), where l is the number of posterior samples. The posterior samples of the c-th regressive coefficient on the d-th landmark for dimension g is `betaOUT[:,(c-1)*k + d, g]`;
 - `sigmaOUT`: dimensions (l,k, k), where l is the number of posterior samples. The element `sigmaOUT[j,:,:]` contains the j-th posterior sample of the covariance matrix;
 - `rmatOUT`: dimensions (l,p, p,n), where l is the number of posterior samples. The element `rmatOUT[j,:,:,a]` contains the j-th posterior sample of the a-th rotation matrix.
 

## Citing

See `CITATION.bib` or use

"Di Noia, A., Mastrantonio, G., and Jona Lasinio, G., “Bayesian Size-and-Shape regression modelling”, <i>arXiv e-prints</i>, 2023".






