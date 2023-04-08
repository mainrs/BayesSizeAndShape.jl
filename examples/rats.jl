using Random
using Distributions
using LinearAlgebra
using StatsBase
using Kronecker 
using DataFrames
using StatsModels 
using CategoricalArrays
using Plots
using BayesSizeAndShape

dataset_rats = dataset("rats");
dataset_desciption("rats")
dataset_names()

landmark = dataset_rats.x;
landmark = landmark ./ 100.0

plot(landmark[:,1,1], landmark[:,2,1],legend = false, color = cgrad(:tab20, 21)[Int64(subject[1])])
for i = 2:size(landmark,3)
    plot!(landmark[:,1,i], landmark[:,2,i], color = cgrad(:tab20, 21)[Int64(subject[i])])
end
title!("Landmarks")


sizeshape = sizeshape_helmertproduct_reflection(landmark);
plot(sizeshape[:,1,1], sizeshape[:,2,1],legend = false, color = cgrad(:tab20, 21)[Int64(subject[1])])
for i = 2:size(landmark,3)
    plot!(sizeshape[:,1,i], sizeshape[:,2,i], color = cgrad(:tab20, 21)[Int64(subject[i])])
end
title!("Size And Shape")

subject = dataset_rats.no;
time = dataset_rats.time;
covariates = DataFrame(
    time = time,
    subject = categorical(string.(subject))
);

outmcmc = SizeAndShapeWithReflectionMCMC(
    landmark,
    @formula(1 ~ 1+time + subject),
    covariates,
    (iter=1000, burnin=200, thin=2),
    Normal(0.0,100000.0),#
    InverseWishart(k + 2, 5.0 * Matrix{Float64}(I, k, k))
);

betaout = posterior_samples_beta(outmcmc);
sigmaout = posterior_samples_sigma(outmcmc);

predictive_mean = sample_predictive_zbr(outmcmc);
predictive_obs = sample_predictive_zbr_plus_epsilon(outmcmc);