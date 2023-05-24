#### #### #### #### #### #### 
#### BETA
#### #### #### #### #### #### 

function sampler_beta(xdata::Array{Float64,3}, designmatrix::Array{Float64,3}, prior::NoPriorBeta, betaMCMC::Matrix{Float64}, sigmaMCMC::Matrix{Float64}, valp::Valuep, betastandMCMC::Matrix{Float64})

end
function sampler_beta(xdata::Array{Float64,3}, designmatrix::Array{Float64,3}, prior::Normal, betaMCMC::Matrix{Float64}, sigmaMCMC::Matrix{Float64}, valp::Valuep, betastandMCMC::Matrix{Float64})

    kd::Int64 = size( betaMCMC,1)
    p::Int64 = size( betaMCMC,2)
    n::Int64 = size(xdata,3)

    invMat = Symmetric(inv(sigmaMCMC))
    for ip = 1:p

        Vp::Matrix{Float64} = zeros(Float64, kd, kd)
        Mp::Vector{Float64} = zeros(Float64, kd)

        Vp[:, :] = Diagonal([1.0 / params(prior)[2]^2 for i = 1:kd])
        Mp[:] .= params(prior)[1] / params(prior)[2]^2

        for j = 1:n
            Vp[:, :] += transpose(designmatrix[:, :, j]) * invMat * designmatrix[:, :, j]
            Mp[:, :] += transpose(designmatrix[:, :, j]) * invMat * xdata[:,ip,j]
        end

        Vp = Symmetric(inv(Vp))
        Mp = Vp*Mp

        betaMCMC[:,ip] = rand(MvNormal(Mp,Vp))
        betastandMCMC[:, ip] = betaMCMC[:, ip]
       
    end
    
    #standardize_reg(betastandMCMC, valp)

end

#### #### #### #### #### #### 
#### SIGMA
#### #### #### #### #### #### 

function sampler_sigma(xdata::Array{Float64,3}, designmatrix::Array{Float64,3}, prior::NoPriorSigma, betaMCMC::Matrix{Float64}, sigmaMCMC::Matrix{Float64}, sigmatype::SigmaType)

end
function sampler_sigma(xdata::Array{Float64,3}, designmatrix::Array{Float64,3}, prior::InverseWishart, betaMCMC::Matrix{Float64}, sigmaMCMC::Matrix{Float64}, sigmatype::GeneralSigma)

    kd::Int64 = size(betaMCMC, 1)
    p::Int64 = size(betaMCMC, 2)
    n::Int64 = size(xdata, 3)

    #invMat = Symmetric(inv(sigmaMCMC))


    nup::Float64 = params(prior)[1] + Float64(n*p) 
    Psip = deepcopy(params(prior)[2].mat)
    for ip = 1:p
        error("rivedere")
        for j = 1:n
            app = xdata[:, p, j] - designmatrix[:,:,j] * betaMCMC[:,p]
            Psip[:,:] += app*transpose(app)
        end

        

        sigmaMCMC[:, :] = rand(InverseWishart(nup, Symmetric(Psip).data))
    end

end



#### #### #### #### #### #### 
#### rmat
#### #### #### #### #### #### 

function sampler_rmat(xdata::Array{Float64,3}, designmatrix::Array{Float64,3}, betaMCMC::Matrix{Float64}, sigmaMCMC::Matrix{Float64}, angleMCMC::Matrix{Float64}, ydata::Array{Float64,3}, valp::Valuep, reflection::KeepReflection, rmatMCMC::Array{Float64,3}, samp_rmat::donotsamplermat)


end
function sampler_rmat(xdata::Array{Float64,3}, designmatrix::Array{Float64,3}, betaMCMC::Matrix{Float64}, sigmaMCMC::Matrix{Float64}, angleMCMC::Matrix{Float64}, ydata::Array{Float64,3}, valp::Valuep2, reflection::KeepReflection, rmatMCMC::Array{Float64,3}, samp_rmat::dosamplermat)


    
    kd::Int64 = size(betaMCMC, 1)
    p::Int64 = size(betaMCMC, 2)
    n::Int64 = size(xdata, 3)
    k::Int64 = size(xdata,1)

    invMat = Symmetric(inv(sigmaMCMC))
    for j = 1:n
    
        mui = designmatrix[:, :, j] * betaMCMC[:, :]
        #@toggled_assert size(mui) == (k,p)
        
        Ai = transpose(mui) * invMat * ydata[:,:,j]
        #println(size(Ai))

        x1 = Ai[1, 1] + Ai[2, 2]
        x2 = Ai[2, 1] - Ai[1, 2]
        #x2 = Ai[1, 2] - Ai[2, 1]

        kvonmises = sqrt(x1^2 + x2^2)

        muvonmises = atan(x2 / kvonmises, x1 / kvonmises)
        #muvonmises = atan(x1 / kvonmises, x2 / kvonmises)

        angleMCMC[1, j] = rand(VonMises(muvonmises, kvonmises))

        rmatMCMC[1, 1, j] = cos(angleMCMC[1, j])
        rmatMCMC[1, 2, j] = -sin(angleMCMC[1, j])
        rmatMCMC[2, 1, j] = sin(angleMCMC[1, j])
        rmatMCMC[2, 2, j] = cos(angleMCMC[1, j])


        #xdata[:, :, j] = ydata[:, :, j] * transpose(rmatMCMC[:,:,j])
    end
    compute_xdata(xdata, ydata, rmatMCMC)

   
   

end

function sampler_rmat(xdata::Array{Float64,3}, designmatrix::Array{Float64,3}, betaMCMC::Matrix{Float64}, sigmaMCMC::Matrix{Float64}, angleMCMC::Matrix{Float64}, ydata::Array{Float64,3}, valp::Valuep3, reflection::KeepReflection, rmatMCMC::Array{Float64,3}, samp_rmat::dosamplermat)

    error("sampler_rmat p=3 - Not completed yet")

    kd::Int64 = size(betaMCMC, 1)
    p::Int64 = size(betaMCMC, 2)
    n::Int64 = size(xdata, 3)
    k::Int64 = size(xdata, 1)

    invMat = Symmetric(inv(sigmaMCMC))
    for j = 1:n

        mui = designmatrix[:, :, j] * betaMCMC[:, :]
        #@toggled_assert size(mui) == (k,p)

        Ai = transpose(mui) * invMat * ydata[:, :, j]
        #println(size(Ai))

        


        x1 = Ai[1, 1] + Ai[2, 2]
        x2 = Ai[2, 1] - Ai[1, 2]
        #x2 = Ai[1, 2] - Ai[2, 1]

        kvonmises = sqrt(x1^2 + x2^2)

        muvonmises = atan(x2 / kvonmises, x1 / kvonmises)
        #muvonmises = atan(x1 / kvonmises, x2 / kvonmises)

        angleMCMC[1, j] = rand(VonMises(muvonmises, kvonmises))

        rmatMCMC[1, 1, j] = cos(angleMCMC[1, j])
        rmatMCMC[1, 2, j] = -sin(angleMCMC[1, j])
        rmatMCMC[2, 1, j] = sin(angleMCMC[1, j])
        rmatMCMC[2, 2, j] = cos(angleMCMC[1, j])


        #xdata[:, :, j] = ydata[:, :, j] * transpose(rmatMCMC[:,:,j])
    end
    compute_xdata(xdata, ydata, rmatMCMC)




end