### ### ### ### ### 
### packages
### ### ### ### ### 

using Pkg
Pkg.activate("/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/gitrepo/BayesSizeAndShape/todelete/sim")
using Random, Distributions, LinearAlgebra, StatsBase
using DataFrames, StatsModels, CategoricalArrays
using Revise
using BayesSizeAndShape
using RCall

### ### ### ### ### 
### simulations
### ### ### ### ### 

n::Int64 = 100;
p::Int64 = 3;
k::Int64 = 10;
d::Int64 = 3;



# regression
reg::Matrix{Float64} = zeros(Float64, k*d, p);
eg::Matrix{Float64} = zeros(Float64, k*d, p);
reg[:] = rand(Normal(0.0, 1.0), prod(size(reg)));
constadd = 10.0
reg[1,:] = [0.0,0.0,0.0].+ constadd
reg[2,:] = [10.0,0.0,-10].+ constadd
reg[3,:] = [20.0,10.0,-10].+ constadd
reg[4,:] = [20.0,20.0,-10].+ constadd
reg[5,:] = [10.0,20.0,-10].+ constadd
reg[6,:] = [0.0,20.0,-20].+ constadd
reg[7,:] = [-10.0,20.0,-20].+ constadd
reg[8,:] = [-20.0,20.0,-20].+ constadd
reg[9,:] = [-20.0,10.0,-20].+ constadd
reg[10,:] = [-10.0,10.0,-20].+ constadd

BayesSizeAndShape.standardize_reg(reg::Matrix{Float64}, BayesSizeAndShape.ValueP3(), BayesSizeAndShape.GramSchmidtMean());
#BayesSizeAndShape.standardize_reg_computegamma(reg::Matrix{Float64}, BayesSizeAndShape.ValueP3(), BayesSizeAndShape.GramSchmidtMean())

zmat = DataFrame(
    x1 = rand(Normal(10.0,1.0 ),n),
    x2 = sample(["A", "B"],n)
)
zmat[:,1] = (zmat[:,1] .- mean(zmat[:,1])) ./ std(zmat[:,1])
zmat.x2 = categorical(zmat.x2)
zmat_modmat_ModelFrame = ModelFrame(@formula(1 ~ 1+x1+x2), zmat);
zmat_modmat = ModelMatrix(zmat_modmat_ModelFrame).m
design_matrix = BayesSizeAndShape.compute_designmatrix(zmat_modmat, k); # dimensions  k, k * d, n


# covariance
sigma::Symmetric{Float64,Matrix{Float64}} = Symmetric(rand(InverseWishart(k + 2, 5.0 * Matrix{Float64}(I, k, k))));


dataset_complete = zeros(Float64,k,p,n);
#dataset = zeros(Float64, k, p, n);
for i_n = 1:n
    for i_p = 1:p
        dataset_complete[:, i_p, i_n] = rand(MvNormal(design_matrix[:, :, i_n] * reg[:, i_p], sigma))
        
    end
end

## TESTS
helmmat = BayesSizeAndShape.RemoveLocationHelmert(k, BayesSizeAndShape.ValueP3());
dataset_complete_landmark = zeros(Float64,k+1,p,n);
for i = 1:n
    dataset_complete_landmark[:,:,i] = transpose(helmmat.matrix)*dataset_complete[:,:,i] .+ 1/(k+1)
end
helmdata = BayesSizeAndShape.remove_location(dataset_complete_landmark, helmmat);


maximum(helmdata[:,:,:]-dataset_complete[:,:,:])
minimum(helmdata[:,:,:]-dataset_complete[:,:,:])

ssdata, ssdata_rotmat = BayesSizeAndShape.compute_sizeshape_fromnolocdata(helmdata, BayesSizeAndShape.KeepReflection(), BayesSizeAndShape.ValueP3());

angle_data = Matrix{Float64}(undef,3,n)
for i = 1:size(ssdata,3)
    BayesSizeAndShape.compute_angle_from_rmat(i,angle_data,  ssdata_rotmat, BayesSizeAndShape.ValueP3(), BayesSizeAndShape.KeepReflection())
end


for i in 1:size(ssdata_rotmat,3)
    println(det(ssdata_rotmat[:,:,i]))
    if det(ssdata_rotmat[:,:,i])<0
        error("")
    end
    println(sum( (helmdata[:,:,i] - ssdata[:,:,i]*transpose(ssdata_rotmat[:,:,i])) ))
end


dataset = BayesSizeAndShape.SSDataType(dataset_complete_landmark,  BayesSizeAndShape.KeepReflection(),helmmat,BayesSizeAndShape.ValueP3(),BayesSizeAndShape.DoNotRemoveSize(), BayesSizeAndShape.GramSchmidtMean());

#dataset.nolocdata[:,:,5] -dataset_complete[:,:,5]
#####



### ### ### ### ### 
### MCMC
### ### ### ### ### 
molt::Int64 = 3
mcmcOUT = generalSizeAndShapeMCMC(;
    landmarks = dataset_complete_landmark,
    fm = @formula(1 ~ 1+ x1 + x2),
    covariates = zmat,
    iterations=(iter=1000*molt, burnin=200*molt, thin=2*molt),

    betaprior = Normal(0.0,100000.0),#
    sigmaprior=InverseWishart(k + 2, 5.0 * Matrix{Float64}(I, k, k)),
    rmatdosample = true,

    #betaprior = BayesSizeAndShape.NoPriorBeta(),
    #sigmaprior = BayesSizeAndShape.NoPriorSigma(),
    #rmatdosample = false,

    #beta_init= zeros(Float64, k * d, p),
    #sigma_init = Symmetric(Matrix{Float64}(I, k, k)),
    #rmat_init = reshape(vcat([Matrix{Float64}(I, p, p)[:] for i = 1:n]...), (p, p, n)),

    beta_init= reg,
    sigma_init = sigma,
    rmat_init = ssdata_rotmat,

    meanmodel = ["linear"][1],
    identification = ["gramschmidt"][1],
    keepreflection = ["no", "yes"][2],
    removelocation = ["no", "helmert"][2],
    removesize = ["no", "norm"][1],
    
    verbose= true
);

#mcmcOUT =BayesSizeAndShape.SizeAndShapeWithReflectionMCMC(
#    dataset_complete_landmark,
#    @formula(1 ~ 1+ x1 + x2),
#    zmat,
#    (iter=1000, burnin=200, thin=2),

#    Normal(0.0,100000.0),#
#    InverseWishart(k + 2, 5.0 * Matrix{Float64}(I, k, k))
#);

predictive_mean = sample_predictive_zbr(mcmcOUT);
predictive_obs = sample_predictive_zbr_plus_epsilon(mcmcOUT);

betaOUT = mcmcOUT.posteriorsamples.beta;
sigmaOUT = mcmcOUT.posteriorsamples.sigma;
rmatOUT = mcmcOUT.posteriorsamples.rmat;
angleOUT = mcmcOUT.posteriorsamples.angle;

@rput angle_data;
@rput predictive_mean;
@rput predictive_obs;
@rput betaOUT;
@rput sigmaOUT;
@rput rmatOUT;
@rput angleOUT;
@rput reg;
@rput sigma;
@rput ssdata_rotmat;
@rput ssdata;
@rput n;
@rput p;
@rput k;


R"""

    DIR = "/Users/gianlucamastrantonio/Dropbox (Politecnico di Torino Staff)/lavori/gitrepo/BayesSizeAndShape/todelete/sim/plot/"

    pdf(paste(DIR, "parametersP3.pdf",sep=""))

    par(mfrow=c(3,3))
    for(i in 1:dim(betaOUT)[2])
    {

        plot(betaOUT[,i], type="l")
        abline(h = c(reg)[i], col=2, lwd=2)

    }
    par(mfrow=c(3,3))
    for(i in 1:dim(sigmaOUT)[2])
    {

        plot(sigmaOUT[,i], type="l")
        abline(h = c(sigma)[i], col=2, lwd=2)

    }
    par(mfrow=c(3,3))
    for(i in 1:dim(rmatOUT)[2])
    {

        plot(rmatOUT[,i], type="l")
        abline(h = c(ssdata_rotmat)[i], col=2, lwd=2)

    }
    par(mfrow=c(3,3))
    for(i in 1:dim(angleOUT)[2])
    {

        plot(angleOUT[,i], type="l", main=colnames(angleOUT)[i])
        abline(h = c(angle_data)[i], col=2, lwd=2)
    }
    dev.off()

    pdf(paste(DIR, "predictive_meanP3.pdf",sep=""))

    meanpred = colMeans(predictive_mean)
    par(mfrow=c(2,2))
    for(i in 1:n)
    {
        plot(matrix(meanpred[(i-1)*(k*p) + 1:(k*p)],ncol=2))
        points(ssdata[,,i], col=2)

        plot(c(ssdata[,,i]),c(matrix(meanpred[(i-1)*(k*p) + 1:(k*p)],ncol=2)))
        abline(a=0,b=1,col=2, lwd=2)
    }

    par(mfrow=c(3,3))
    for(i in 1:dim(predictive_mean)[2])
    {

        plot(predictive_mean[,i], type="l", main= c(ssdata)[i])
        abline(h = c(ssdata)[i], col=2, lwd=2)

    }
    
    dev.off()

    pdf(paste(DIR, "predictive_obsP3.pdf",sep=""))

    meanpred = colMeans(predictive_obs)
    par(mfrow=c(2,2))
    for(i in 1:n)
    {
        plot(matrix(meanpred[(i-1)*(k*p) + 1:(k*p)],ncol=2))
        points(ssdata[,,i], col=2)

        plot(c(ssdata[,,i]),c(matrix(meanpred[(i-1)*(k*p) + 1:(k*p)],ncol=2)))
        abline(a=0,b=1,col=2, lwd=2)
    }

    par(mfrow=c(3,3))
    for(i in 1:dim(predictive_obs)[2])
    {

        plot(predictive_obs[,i], type="l", main= c(ssdata)[i])
        abline(h = c(ssdata)[i], col=2, lwd=2)

    }
    
    
    dev.off()

"""


##### TESTS
#if true == true

#    design_matrix
#    maximum(abs.(mcmcOUT.meantype.designmatrix-design_matrix))


#    zmat_modmat - mcmcOUT.meantype.model_matrix
#end








