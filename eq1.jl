function eq1(N::Int, J::Int, pi_nij::Array{Float64, 3}, hat_kappa_nij::Array{Float64, 3}, hat_x_ij::Array{Float64, 2}, theta_j::Array{Float64, 1}, hat_T_nj::Array{Float64, 2}, gamma_nj::Array{Float64, 2})
    a = zeros(N, J)
    for n in 1:N
        for j in 1:J
            a[n, j] = sum(pi_nij[n, :, j] .* (hat_kappa_nij[n, :, j] .* hat_x_ij[:, j]).^(-theta_j[j]) .* (hat_T_nj[:, j]).^(gamma_nj[n, j] .* theta_j[j])).^(-1/theta_j[j])
        end
    end
    return a
end
