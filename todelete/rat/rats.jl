#################
# Packages     ##
#################

using Pkg
Pkg.activate("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/gitrepo/BayesSizeAndShape/todelete/rat")

using Revise
using Random, Distributions, LinearAlgebra,  StatsBase
using Kronecker, DataFrames,StatsModels, CategoricalArrays
using Plots

using BayesSizeAndShape

dataset_rats = dataset("rats");
dataset_desciption("rats")

landmark = dataset_rats.x;
landmark = landmark ./ 100.0
subject = dataset_rats.no;

time = dataset_rats.time;


## design matrix

X = DataFrame(
    time = time,
    subject = categorical(string.(subject))
);

## size and shape data
sizeshape = sizeshape_helmertproduct_reflection(landmark);

k = size(sizeshape,1); 
n = size(sizeshape,3); 
p = size(sizeshape,2); 

## plots of the data



plot(landmark[:,1,1], landmark[:,2,1],legend = false, color = cgrad(:tab20, 21)[Int64(subject[1])])
for i = 2:size(landmark,3)
    plot!(landmark[:,1,i], landmark[:,2,i], color = cgrad(:tab20, 21)[Int64(subject[i])])
end
title!("Landmarks")
plot(sizeshape[:,1,1], sizeshape[:,2,1],legend = false, color = cgrad(:tab20, 21)[Int64(subject[1])])
for i = 2:size(landmark,3)
    plot!(sizeshape[:,1,i], sizeshape[:,2,i], color = cgrad(:tab20, 21)[Int64(subject[i])])
end
title!("Size And Shape")




## covariates


### ### ### ### ### 
### MCMC
### ### ### ### ### 

outmcmc = SizeAndShapeWithReflectionMCMC(
    landmark,
    @formula(1 ~ 1+time + subject),
    X,
    (iter=1000, burnin=200, thin=2),
    Normal(0.0,100000.0),#
    InverseWishart(k + 2, 5.0 * Matrix{Float64}(I, k, k))
);

predictive_mean = sample_predictive_zbr(outmcmc);
predictive_obs = sample_predictive_zbr_plus_epsilon(outmcmc);

betaOUT = outmcmc.posteriorsamples.beta;
sigmaOUT = outmcmc.posteriorsamples.sigma;
rmatOUT = outmcmc.posteriorsamples.rmat;

describe(betaOUT[:,1:10], :all)

#@rput betaOUT;
#@rput sigmaOUT;



#R" save(betaOUT, sigmaOUT , file='provarats.RData')"

#if flag == 1

#    prediction = BayesSizeAndShape.predictmean(outmcmc);
    
#    @rput dataset;
#    @rput prediction;

#    R"""
#    #require(emdbook)
#    #require(coda)
#        OUTpred = list()
#        # ogni elemento della lista prediction contiene i posterior sample di un'osservazione
#        for(select_obs in 1:length(prediction))
#        {
            
#            land.m <- apply(prediction[[select_obs]],c(2,3),mean)
#            qq1 <- apply(prediction[[select_obs]][,,1],c(2),quantile,prob=c(0.025,0.975))
#            qq2 <- apply(prediction[[select_obs]][,,2],c(2),quantile,prob=c(0.025,0.975))

#            OUTpred[[select_obs]]<-list(meanland=land.m,CIland1=qq1, CIland2=qq2)
#        }
#        #aggiungere bb<-mcmc(prediction[[select_obs]])
        
#        pdf("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/SizeAndShape/JuliaStuff/GIAN/chainsGIAN.pdf")
#        #theta = 0
#        #R = matrix(c(cos(theta),sin(theta),-sin(theta), cos(theta)), ncol=2)
#        plot(dataset[,,1], xlim=c(-300,1000)/100,ylim=c(-500,500)/100)
#        points(OUTpred[[1]]$meanland, col=2)
#        for(select_obs in 1:length(prediction))
#        {
#            points(dataset[,,select_obs])
#            points(OUTpred[[select_obs]]$meanland, col=2)
#            for(ll in 1:7){
#                bb<-mcmc(prediction[[select_obs]][,ll,])
#                HPDregionplot(bb,add=T,col=2,lty=2)
#                }
#        }
        
        
#        for(p in 1:dim(betaOUT)[3])
#        {
#            par(mfrow=c(3,3))
#            for(i in 1:dim(betaOUT)[2])
#            {
#                plot(betaOUT[,i,p], type="l")
#            } 
#        }
#        dev.off()
#        #}
#        #OUTpred[[select_obs]]<-list(meanland=land.m,CIland=qq.m)
#    """
#    # mean_pred contiene i campioni a posteriori
#    R" save(OUTpred, file='OUTpred1.RData')"
#    else print("Done posterior estimation")
#end