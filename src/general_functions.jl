

function compute_designmatrix(zmat::Matrix{Float64}, k::Int64)::Array{Float64,3}

    n, d = size(zmat)
    res::Array{Float64,3} = zeros(Float64, k, k * d, n)
    for i = 1:n
        res[:, :, i] = kronecker(transpose(zmat[i, :]), Matrix{Float64}(I, k, k))
    end

    return res

end

function compute_designmatrix(zmat::DataFrame, k::Int64)::Array{Float64,3}

    n, d = size(zmat)
    res::Array{Float64,3} = zeros(Float64, k, k * d, n)
    for i = 1:n
        res[:, :, i] = kronecker(transpose([values(zmat[i, :])...]), Matrix{Float64}(I, k, k))
    end

    return res

end

function compute_xdata(xdata::Array{Float64,3}, ydata::Array{Float64,3}, rmat::Array{Float64,3})

    for i = axes(xdata,3)
        xdata[:, :, i] = ydata[:, :, i] * transpose(rmat[:, :, i])
    end

    return nothing

end


function compute_ss_from_pre_p3_keep_reflection(xdata::Array{Float64,3}, ydata::Array{Float64,3}, ret::Array{Float64,3})

    for i = axes(xdata, 3)

        app = svd(xdata[:, :, i])
        U = app.U
        V = transpose(app.Vt)
        if det(V) > 0

            ydata[:, :, i] = U * Diagonal(app.S)
            ret[:, :, i] = V


        else



            U[:, 3] = -1.0 * U[:, 3]
            V[:, 3] = -1.0 * V[:, 3]


            ydata[:, :, i] = U * Diagonal(app.S)
            ret[:, :, i] = V


        end

    end

end

function compute_ss_from_pre_p3_nokeep_reflection(xdata::Array{Float64,3}, ydata::Array{Float64,3}, ret::Array{Float64,3})

    for i = axes(xdata, 3)

        app = svd(xdata[:, :, i])

        ydata[:, :, i] = app.U * Diagonal(app.S)
        ret[:, :, i] = transpose(app.Vt)

    end

end



function compute_ss_from_pre_p2_keep_reflection(xdata::Array{Float64,3}, ydata::Array{Float64,3}, ret::Array{Float64,3})

    for i = axes(xdata, 3)
    
        app = svd(xdata[:, :, i])
        U = app.U
        V = transpose(app.Vt)
        if sign(V[1, 2]) != sign(V[2, 1])

            ydata[:, :, i] = U * Diagonal(app.S)
            ret[:, :, i] = V


        else

            
            
            U[:,2] = -1.0*U[:,2]
            V[:, 2] = -1.0 * V[:, 2]


            ydata[:, :, i] = U * Diagonal(app.S)
            ret[:, :, i] = V
            

        end

    end

end

function compute_ss_from_pre_p2_nokeep_reflection(xdata::Array{Float64,3}, ydata::Array{Float64,3}, ret::Array{Float64,3})

    for i = axes(xdata, 3)

        app = svd(xdata[:, :, i])

        ydata[:, :, i] = app.U * Diagonal(app.S)
        ret[:, :, i] = transpose(app.Vt)

    end

end


"""
    compute_ss_from_pre(xdata::Array{Float64,3},ydata::Array{Float64,3}, keep_reflection::Bool)::Array{Float64,3}

Given the array of preform matrix xdata (xdata[:,:,i] is the i-th preform matrix), the function computes the size-and-shape data and store it in ydata.
The function return the array of associated rotation matrix R, such that xdata = ydata*R. keep_reflection is a boolean that is used to indicate if R in O(p) (keep_reflection=false) or R in SO(p) (keep_reflection=true) 

    compute_ss_from_pre(xdata::Array{Float64,3}, keep_reflection::Bool)::Array{Float64,3}

Given the array of preform matrix xdata (xdata[:,:,i] is the i-th preform matrix), the function computes the size-and-shape data.
The function retunn the size-and-shape ydata. keep_reflection is a boolean that is used to indicate if R in O(p) (keep_reflection=false) or R in SO(p) (keep_reflection=true) 

"""
function compute_ss_from_pre(xdata::Array{Float64,3},ydata::Array{Float64,3}, keep_reflection::Bool)::Array{Float64,3}

    p = size(xdata, 2) 
    ret = zeros(Float64, p, p, size(xdata, 3))
    if (p == 2) & (keep_reflection == true)
        compute_ss_from_pre_p2_keep_reflection(xdata, ydata, ret)
    end

    if (p == 2) & (keep_reflection == false)
        compute_ss_from_pre_p2_nokeep_reflection(xdata, ydata, ret)
    end

    if (p == 3) & (keep_reflection == true)
        compute_ss_from_pre_p3_keep_reflection(xdata, ydata, ret)
    end
    if (p == 3) & (keep_reflection == false)
        compute_ss_from_pre_p3_nokeep_reflection(xdata, ydata, ret)
    end

    if p >3
        error("size(xdata, 2) must be 2 or 3")
    end
    
    return ret

end

function compute_ss_from_pre(xdata::Array{Float64,3}, keep_reflection::Bool)::Array{Float64,3}

    ydata::Array{Float64,3} = deepcopy(xdata)
    p = size(xdata, 2) 
    ret = zeros(Float64, p, p, size(xdata, 3))
    if (p == 2) & (keep_reflection == true)
        compute_ss_from_pre_p2_keep_reflection(xdata, ydata, ret)
    end

    if (p == 2) & (keep_reflection == false)
        compute_ss_from_pre_p2_nokeep_reflection(xdata, ydata, ret)
    end

    if (p == 3) & (keep_reflection == true)
        compute_ss_from_pre_p3_keep_reflection(xdata, ydata, ret)
    end
    if (p == 3) & (keep_reflection == false)
        compute_ss_from_pre_p3_nokeep_reflection(xdata, ydata, ret)
    end
    
    return ydata

end


function compute_angle_from_rmat(i::Int64,angle::Matrix{Float64},  rmat::Array{Float64,3}, valp::Valuep2, reflection::KeepReflection)

    angle[1, i] = atan(rmat[2, 1, i], rmat[1, 1, i])

    return nothing
end

function compute_angle_from_rmat(i::Int64, angle::Matrix{Float64}, rmat::Array{Float64,3}, valp::Valuep3, reflection::KeepReflection)
    #x convention
    angle[1, i] = atan(rmat[3, 1, i], rmat[3, 2, i])
    angle[2, i] = acos(rmat[3, 3, i])
    angle[3, i] = -atan(rmat[1, 3, i], rmat[2, 3, i])

    return nothing
end


function standardize_reg(reg::Matrix{Float64}, valp::Valuep2)

    a1 = reg[1, :]
    a2 = reg[2, :]
    g1 = a1 / norm(a1, 2)
    g2 = (a2 - (transpose(a2) * g1) * g1) / norm(a2 - (transpose(a2) * g1) * g1, 2)
    gammamat = reshape([g1; g2], (2, 2))

    if det(gammamat) < 0
        gammamat[:, 2] = - gammamat[:, 2]
    end

    reg[:,:] = reg * gammamat

end

function standardize_reg(reg::Matrix{Float64}, valp::Valuep3)

    a1 = reg[1, :]
    a2 = reg[2, :]
    a3 = reg[3, :]
    g1 = a1 / norm(a1, 2)
    g2 = (a2 - (transpose(a2) * g1) * g1) / norm(a2 - (transpose(a2) * g1) * g1, 2)
    g3 = (a3 - (transpose(a3) * g1) * g1 - (transpose(a3) * g2) * g2) / norm(a3 - (transpose(a3) * g1) * g1 - (transpose(a3) * g2) * g2, 2)



    gammamat = reshape([g1; g2; g3], (3, 3))

    if det(gammamat) < 0
        gammamat[:, 3] = -gammamat[:, 3]
    end

    reg[:, :] = reg * gammamat

end

"""
    compute_helmertized_configuration(landmark::Array{Float64,3})::Array{Float64,3}

Given the array of landmarks (landmark[:,:,i] is the matrix of landmarks of the i-th shape), it computes and returns the helmertized configuration
"""
function compute_helmertized_configuration(landmark::Array{Float64,3})::Array{Float64,3} 

    k::Int64 = size(landmark,1)
    H::Matrix{Float64} = zeros(Float64,k,k)
    ret::Array{Float64,3} = zeros(Float64,size(landmark,1)-1,size(landmark,2),size(landmark,3))
    for i = 1:k
        H[i,i] = (i-1)/(i*(i-1))
    end
    for i = 2:k
        for j = 1:(i-1)
            H[i,i] = -1.0/(i*(i-1))
        end
    end
    H = H[2:end,:]

    for i = 1:size(landmark,3)
        ret[:,:,i] = H*landmark[:,:,i]
    end

    return ret

end