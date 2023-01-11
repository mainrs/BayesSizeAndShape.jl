
function mcmc(;
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
    rmat_init::Array{Float64,3} ,
    dormat::Bool,
    reflection::Reflection = KeepReflection(),
    sigmatype::SigmaType = GeneralSigma()
)


    
    #### #### #### #### #### 
    #### check dimensions  and base objects
    #### #### #### #### #### 

    k::Int64 = size(dataset, 1)
    p::Int64 = size(dataset, 2)
    n::Int64 = size(dataset, 3)

    valp = Valuep2()
    if p ==2
        valp = Valuep2()
    end
    if p == 3
        valp = Valuep3()
    end

    designmatrix_v2 = modelmatrix(fm, covariates) #n, d

    @assert sum(designmatrix_v2[:, 1] .== 1) == n "intercept needed"
    for i = 2:size(designmatrix_v2, 2)

        designmatrix_v2[:, i] = designmatrix_v2[:, i] .- mean(designmatrix_v2[:, i])

    end

    designmatrix = compute_dessignmatrix(modelmatrix(fm, covariates), k) # dimensions  k, k * d, n
    
    d::Int64 = size(designmatrix_v2, 2)
    
    @assert size(dataset, 3) == size(covariates,1) # n
    @assert size(dataset, 3) == size(designmatrix_v2,1) # n
    @assert size(designmatrix, 1) == size(dataset, 1) # k
    @assert size(designmatrix, 3) == size(dataset, 3) # n
    @assert size(designmatrix, 2) == size(dataset, 1) * size(designmatrix_v2,2) #k d
    @assert size(beta_init, 2) == size(dataset, 2) # p
    @assert size(beta_init, 1) == size(designmatrix, 2) # k d
    @assert size(sigma_init, 1) == size(dataset, 1)
    @assert size(rmat_init, 1) == size(dataset, 2)
    @assert size(rmat_init, 2) == size(dataset, 2)
    @assert size(rmat_init, 3) == size(dataset, 3)

    designmatrix_v2 = nothing

    xdata = deepcopy(dataset)
    compute_xdata(xdata, dataset, rmat_init)

    

    #### #### #### #### #### 
    #### MCMC object
    #### #### #### #### #### 

    betaMCMC::Matrix{Float64} = deepcopy(beta_init)
    betastandMCMC::Matrix{Float64} = deepcopy(beta_init)
    sigmaMCMC::Matrix{Float64} = deepcopy(sigma_init)
    rmatMCMC::Array{Float64,3} = deepcopy(rmat_init)

    angleMCMC::Matrix{Float64} = zeros(Float64, ifelse(p==2,1,3),n)

    if typeof(reflection) == KeepReflection
        
        if p == 2
            for i = 1:n
                compute_angle_from_rmat(i, angleMCMC, rmatMCMC, valp, reflection)
            end
        else
            error()
        end

    elseif typeof(reflection) == donotKeepReflection
        error()
    end

    samp_rmat = dosamplermat()
    if dormat == false

        samp_rmat = donotsamplermat()
        
    end
    #### #### #### #### #### 
    #### MCMC
    #### #### #### #### #### 
    iter = Int64(iterations.iter)
    burnin = Int64(iterations.burnin)
    thin = Int64(iterations.thin)
    sampletosave = trunc(Int64, round((iter - burnin) / thin))

    #### #### #### #### #### 
    #### MCMC out
    #### #### #### #### #### 

    betaOUT::Array{Float64,3} = zeros(Float64, Tuple([sampletosave; size(betaMCMC)...]))
    sigmaOUT::Array{Float64,3} = zeros(Float64, Tuple([sampletosave; size(sigmaMCMC)...]))
    rmatOUT::Array{Float64,4} = zeros(Float64, Tuple([sampletosave; size(rmatMCMC)...]))
    angleOUT::Array{Float64,3} = zeros(Float64, Tuple([sampletosave; size(angleMCMC)...]))
    # sigmaMCMC::Matrix{Float64} = deepcopy(sigma_init)
    # rmatMCMC::Array{Float64,3} = deepcopy(rmat_init)

    # angleMCMC::Matrix{Float64} = zeros(Float64, ifelse(p == 2, 1, 3), n)



    #### #### #### #### #### 
    #### algrothm
    #### #### #### #### #### 
    iterMCMC = Int64(0)
    thinburnin = burnin
    p1 = Progress(burnin, desc = "burnin ", offset = 0, showspeed = true)
    p2 = Progress(
        burnin + (sampletosave - 1) * thin,
        desc = "iterations ",
        offset = 0,
        showspeed = true,
    )
    isburn = true

    println("Iterations: ", iter)
    println("burnin: ", burnin)
    println("thin: ", thin)
    println("number of posterior samples: ", sampletosave)
    println("Number of threads: ", Threads.nthreads())
    for iMCMC = 1:sampletosave

        for jMCMC = 1:thinburnin

            iterMCMC += 1
            ProgressMeter.next!(p2; showvalues = [(:iterations, iterMCMC)])



            sampler_beta(xdata, designmatrix, betaprior, betaMCMC, sigmaMCMC, valp, betastandMCMC)
            betaMCMC[:, :] = betastandMCMC[:,:]
            sampler_rmat(xdata, designmatrix, betaMCMC, sigmaMCMC, angleMCMC, dataset, valp, reflection, rmatMCMC, samp_rmat)
            sampler_sigma(xdata, designmatrix, sigmaprior, betaMCMC, sigmaMCMC, sigmatype)
            


        end

        thinburnin = thin
        isburn = false

        
        betaOUT[iMCMC, :, :] = betastandMCMC[:, :]
        sigmaOUT[iMCMC, :, :] = sigmaMCMC[:, :]
        rmatOUT[iMCMC, :, :, :] = rmatMCMC[:, :, :]
        angleOUT[iMCMC, :, :] = angleMCMC[:, :]

    end

    return betaOUT, sigmaOUT, rmatOUT, angleOUT


end
