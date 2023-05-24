
function sim(idata::Int64, NAMEGEN = "SIM",ndata::Int64, nsim::Int64,betainsideCI::Matrix{Float64} , sigmainsideCI::Matrix{Float64},betalengthCI::Matrix{Float64}, sigmalengthCI::Matrix{Float64})

    for ip = 1:2
        for isim = 1:nsim
            
    
            n::Int64 = [20,100,20,100, 20,100,20,100][isim];
            p::Int64 = [2,3][ip];
            k::Int64 = [10,10,20,20, 10,10,20,20][isim];
            d::Int64 = 3;
            var::Float64 = [1.0,1.0,1.0,1.0, 10.0,10.0,10.0,10.0][isim]
    
            NAME = NAMEGEN * "_" * string(n)* "_" * string(p) * "_" * string(k) * "_" * string(var)
    
            
    
    
            # regression
            reg::Matrix{Float64} = zeros(Float64, k*d, p);
            #eg::Matrix{Float64} = zeros(Float64, k*d, p);
            reg[:] = rand(Normal(5.0, 1.0), prod(size(reg)));
            #constadd = 10.0
            #if p == 2
            #    reg[1,:] = [0.0,0.0].+ constadd
            #    reg[2,:] = [10.0,0.0].+ constadd
            #    reg[3,:] = [20.0,10.0].+ constadd
            #    reg[4,:] = [20.0,20.0].+ constadd
            #    reg[5,:] = [10.0,20.0].+ constadd
            #    reg[6,:] = [0.0,20.0].+ constadd
            #    reg[7,:] = [-10.0,20.0].+ constadd
            #    reg[8,:] = [-20.0,20.0].+ constadd
            #    reg[9,:] = [-20.0,10.0].+ constadd
            #    reg[10,:] = [-10.0,10.0].+ constadd
            #else
            #    reg[1,:] = [0.0,0.0,0.0].+ constadd
            #    reg[2,:] = [10.0,0.0,-10].+ constadd
            #    reg[3,:] = [20.0,10.0,-10].+ constadd
            #    reg[4,:] = [20.0,20.0,-10].+ constadd
            #    reg[5,:] = [10.0,20.0,-10].+ constadd
            #    reg[6,:] = [0.0,20.0,-20].+ constadd
            #    reg[7,:] = [-10.0,20.0,-20].+ constadd
            #    reg[8,:] = [-20.0,20.0,-20].+ constadd
            #    reg[9,:] = [-20.0,10.0,-20].+ constadd
            #    reg[10,:] = [-10.0,10.0,-20].+ constadd
            #end
    
    
            if p == 2
                BayesSizeAndShape.standardize_reg(reg::Matrix{Float64}, BayesSizeAndShape.ValueP2(), BayesSizeAndShape.GramSchmidtMean());
            else
                BayesSizeAndShape.standardize_reg(reg::Matrix{Float64}, BayesSizeAndShape.ValueP3(), BayesSizeAndShape.GramSchmidtMean());
            end
    
    
    
            zmat = DataFrame(
                x1 = rand(Normal(10.0,1.0 ),n),
                x2 = sample(["A", "B"],n)
            )
            zmat[:,1] = (zmat[:,1] .- mean(zmat[:,1])) ./ std(zmat[:,1])
            zmat.x2 = categorical(zmat.x2)
            zmat_modmat_ModelFrame = ModelFrame(@formula(1 ~ 1+x1 + x2 ), zmat);
            zmat_modmat = ModelMatrix(zmat_modmat_ModelFrame).m
            design_matrix = BayesSizeAndShape.compute_designmatrix(zmat_modmat, k); # dimensions  k, k * d, n
    
    
            # covariance
            sigma::Symmetric{Float64,Matrix{Float64}} = Symmetric(rand(InverseWishart(k + 2, 5.0 * Matrix{Float64}(I, k, k))));
            sigma.data[:,:] = sigma.data[:,:] .* var 
    
            dataset_complete = zeros(Float64,k,p,n);
            #dataset = zeros(Float64, k, p, n);
            for i_n = 1:n
                for i_p = 1:p
                    dataset_complete[:, i_p, i_n] = rand(MvNormal(design_matrix[:, :, i_n] * reg[:, i_p], sigma))
                    
                end
            end
    
    
            if p == 2
                helmmat = BayesSizeAndShape.RemoveLocationHelmert(k, BayesSizeAndShape.ValueP2());
            else
                helmmat = BayesSizeAndShape.RemoveLocationHelmert(k, BayesSizeAndShape.ValueP3());
            end
    
            dataset_complete_landmark = zeros(Float64,k+1,p,n);
            for i = 1:n
                dataset_complete_landmark[:,:,i] = transpose(helmmat.matrix)*dataset_complete[:,:,i] .+ 1/(k+1)
            end
            helmdata = BayesSizeAndShape.remove_location(dataset_complete_landmark, helmmat);
    
    
            if p == 2
                ssdata, ssdata_rotmat = BayesSizeAndShape.compute_sizeshape_fromnolocdata(helmdata, BayesSizeAndShape.KeepReflection(), BayesSizeAndShape.ValueP2());
            else
                ssdata, ssdata_rotmat = BayesSizeAndShape.compute_sizeshape_fromnolocdata(helmdata, BayesSizeAndShape.KeepReflection(), BayesSizeAndShape.ValueP3());
            end
    
    
            if p == 2
                dataset = BayesSizeAndShape.SSDataType(dataset_complete_landmark,  BayesSizeAndShape.KeepReflection(),helmmat,BayesSizeAndShape.ValueP2(),BayesSizeAndShape.DoNotRemoveSize(), BayesSizeAndShape.GramSchmidtMean());
            else
                dataset = BayesSizeAndShape.SSDataType(dataset_complete_landmark,  BayesSizeAndShape.KeepReflection(),helmmat,BayesSizeAndShape.ValueP3(),BayesSizeAndShape.DoNotRemoveSize(), BayesSizeAndShape.GramSchmidtMean());
            end
    
    
            #dataset.nolocdata[:,:,5] -dataset_complete[:,:,5]
            #####
    
    
    
            ### ### ### ### ### 
            ### MCMC
            ### ### ### ### ### 
            molt::Int64 = 30
            outmcmc = SizeAndShapeWithReflectionMCMC(
                dataset_complete_landmark,
                @formula(1 ~ 1+ x1 + x2 ),
                zmat,
                (iter=3000*molt, burnin=1000*molt, thin=1*molt),
                Normal(0.0,100000.0),#
                InverseWishart(k + 2, 2.0 * Matrix{Float64}(I, k, k))
            );
    
            betaout = posterior_samples_beta(outmcmc);
            sigmaout = posterior_samples_sigma(outmcmc);
    
    
            windex = isim + (ip-1)*nsim
    
            betainsideCI[windex,idata] = 0.0
            npar = 0
            for idim = 1:p
                initpar = [1,2,3][idim]
    
                for ipar in initpar:( d*k)
                    npar += 1
                            
                    ww = ipar + (idim-1)*d*k
                    qq = quantile(betaout[:,ww], [0.025,0.975])
                    if (reg[ipar,idim]>=qq[1]) & (reg[ipar,idim]<=qq[2])
                        
                        betainsideCI[windex,idata] += 1.0
    
                    end
                    betalengthCI[windex,idata] += qq[2]-qq[1]
                    
                end
                    
            end
            betainsideCI[windex,idata] = betainsideCI[windex,idata] / npar
            betalengthCI[windex,idata] = betalengthCI[windex,idata] / npar
    
            sigmainsideCI[windex,idata] = 0.0
            npar = 0
            for ipar = 1:k
    
                for jpar = 1:ipar
    
                    npar += 1
    
                    ww = jpar + (ipar-1)*k
                    qq = quantile(sigmaout[:,ww], [0.025,0.975])
                    
                    if (sigma[ipar,jpar]>=qq[1]) & (sigma[ipar,jpar]<=qq[2])
                        
                        sigmainsideCI[windex,idata] += 1.0
    
                    end
                    sigmalengthCI[windex,idata] += qq[2]-qq[1]
                end
    
            end
            sigmainsideCI[windex,idata] = sigmainsideCI[windex,idata] / npar
            sigmalengthCI[windex,idata] = sigmalengthCI[windex,idata] / npar
    
    
            println(trunc.(betainsideCI[1:windex,idata],digits=2))
            println(trunc.(sigmainsideCI[1:windex,idata],digits=2))
    
            println(trunc.(betalengthCI[1:windex,idata],digits=2))
            println(trunc.(sigmalengthCI[1:windex,idata],digits=2))
    
            @rput betainsideCI
            @rput sigmainsideCI
            @rput betalengthCI
            @rput sigmalengthCI
            #@rput betaout;
            #@rput sigmaout;
            #@rput reg;
            #@rput sigma;
            #@rput ssdata_rotmat;
            #@rput ssdata;
            #@rput n;
            #@rput p;
            #@rput k;
            #@rput NAME;
    
            
    
    
        end
    end


end