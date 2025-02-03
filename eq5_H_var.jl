function eq5_H_var(N, J, alpha_nj, hat_omega_n,hat_L_n,hat_H_n, beta_n, In,Ln, Sn, Sn_l,gamma_njk, pi_nij)
    coefs_A = zeros(N, N, J, J)
    for n in 1:N,j in 1:J,k in 1:J,i in 1:N
                    coefs_A[n, i, k, j] = gamma_njk[n, k, j]*pi_nij[i, n, k]
    end
    # Create the A matrix
    A = I(N*J) - reshape(permutedims(coefs_A, (1, 3, 2, 4)), (N * J, N * J))
    A[isinf.(A)].=0
    A[isnan.(A)].=0
    # Create the c vector
    c = zeros(N, J)
    for n in 1:N
        for j in 1:J
            c[n, j] = alpha_nj[n, j] * (hat_omega_n[n] * hat_H_n[n]^(beta_n[n]) * hat_L_n[n]^(1 - beta_n[n]) * (In[n] * Ln[n] + Sn[n] - Sn_l[n]))
        end
    end
    c = reshape(c, (N * J, 1))
    # Solving the system:
    return reshape(A\c,(N,J))
end