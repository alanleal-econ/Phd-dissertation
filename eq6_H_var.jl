function eq6_H_var(N, J, gamma_nj, pi_nij, X_nj, hat_L_n, hat_H_n,beta_n, In, Ln, Sn)
    num=zeros(N)
    den=zeros(N)
    a=zeros(N)
    for n in 1:N
        num[n] = sum(gamma_nj[n,j]*sum(pi_nij[i,n,j]*X_nj[i,j] for i=1:N) for j=1:J)
        den[n] = hat_H_n[n]^(beta_n[n])*hat_L_n[n]^(1 - beta_n[n]) * (In[n] * Ln[n] + Sn[n])
        a[n]=num[n]/den[n]
    end
    return a
end