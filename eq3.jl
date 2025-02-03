function eq3(N::Int, J::Int, pi_nij::Array{Float64,3}, hat_P_nj::Matrix{Float64},
    hat_kappa_nij::Array{Float64,3}, hat_x_ij::Matrix{Float64},
    theta_j::Vector{Float64}, hat_T_nj::Matrix{Float64},
    gamma_nj::Matrix{Float64})
    a = zeros(N, N, J)
    for n in 1:N, i in 1:N, j in 1:J
        a[n, i, j] = pi_nij[n, i, j] * (hat_P_nj[n, j] / (hat_kappa_nij[n, i, j] * hat_x_ij[i, j]))^theta_j[j] * hat_T_nj[n, j]^(gamma_nj[i, j] * theta_j[j])
    end
    return a
end
