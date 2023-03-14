function compute_dessignmatrix(zmat::Matrix{Float64}, k::Int64)::Array{Float64,3}
    println("compute_dessignmatrix is deprecated. Please use compute_designmatrix instead")
    return compute_designmatrix(zmat, k)

end

function compute_dessignmatrix(zmat::DataFrame, k::Int64)::Array{Float64,3}
    println("compute_dessignmatrix is deprecated. Please use compute_designmatrix instead")
    return compute_designmatrix(zmat, k)

end