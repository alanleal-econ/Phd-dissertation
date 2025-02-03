function eq2(N::Int, J::Int, hat_omega_n::Vector{Float64}, gamma_nj::Matrix{T}, hat_P_nj::Matrix{T}, gamma_njk::Array{T,3}) where T<:AbstractFloat
    a = zeros(T, N, J)
    for n in 1:N
        for j in 1:J
            a[n, j] = hat_omega_n[n]^gamma_nj[n, j] * prod(hat_P_nj[n, :] .^ gamma_njk[n, j, :])
        end
    end
    return a
end
