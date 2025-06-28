# Main code
# Created by Alan Leal. E-mail: prof@alanleal-econ.com
N=70
J=26
hat_kappa_nij=ones(N,N,J)
# EU countries code
br_codes=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27]
eu_codes=[29,30,31,36,37,38,39,40,41,42,43,45,46,47,50,51,54,55,56,58,59,61,62,63,65,66,67]
# Baseline:
hat_kappa_nij[:,:,:].=1
# Run one scenario at a time
# Scenario i - Economy-wide trade cost reduction
hat_kappa_nij[br_codes,eu_codes,1:15].=0.7646
hat_kappa_nij[eu_codes,br_codes,1:15].=0.7646
# Scenario ii - Primary sectors surcharge, with trade cost reduction in all other sectors
hat_kappa_nij[br_codes,eu_codes,4:26].=0.7646
hat_kappa_nij[eu_codes,br_codes,4:26].=0.7646
hat_kappa_nij[br_codes,eu_codes,1:3].=1.5 #50% increase in the trade cost
hat_kappa_nij[eu_codes,br_codes,1:3].=1.5 #50% increase in the trade cost


# Endogenous variables
hat_T_nj=ones(N,J)
hat_S_n=zeros(N,1)
Sn_l=zeros(N,1)

# Loading the required package
using LinearAlgebra
using Missings
using Optim
using NLsolve
using NaNStatistics
using JuMP, Ipopt


# Loading functions:
include("eq1.jl")
include("eq2.jl")
include("eq3.jl")
include("eq4_h_endog_v2.jl")
include("eq5_H_var.jl")
include("eq6_H_var.jl")
# Creating empty objects of interest:
hat_P_nj=ones(N,J)
hat_x_ij=ones(N,J)
pi_nij_l=zeros(N,J)
hat_L_n=ones(N)
X_nk_l=zeros(N,J)
hat_Ub=ones(N-26)
hat_omega_n0=ones(N)
hat_H_n=ones(N,1)
hat_L_n=ones(N,1)
epsilon=exp(100)
cont=1
while epsilon>10^(-3) && cont<=1000
    epsilon1=exp(100)
    while epsilon1>(10^(-7))
        hat_P_nj1=hat_P_nj
        hat_P_nj=eq1(N,J,pi_nij,hat_kappa_nij,hat_x_ij,theta_j,hat_T_nj,gamma_nj)
        hat_P_nj[isinf.(hat_P_nj)].=1;hat_P_nj[isnan.(hat_P_nj)].=1;hat_P_nj[iszero.(hat_P_nj)].=1
        hat_x_ij=eq2(N,J,hat_omega_n0,gamma_nj,hat_P_nj,gamma_njk)
        hat_x_ij[isinf.(hat_x_ij)].=1;hat_x_ij[isnan.(hat_x_ij)].=1;hat_x_ij[iszero.(hat_x_ij)].=1
        epsilon1=nansum((hat_P_nj1-hat_P_nj).^2)
    end

    # Eq. 3 now:
    pi_nij_l=eq3(N,J,pi_nij,hat_P_nj,hat_kappa_nij,hat_x_ij,theta_j,hat_T_nj,gamma_nj)
    pi_nij_l[isinf.(pi_nij_l)].=1;pi_nij_l[isnan.(pi_nij_l)].=1
    # Eq. 4:
    # Creating object of interest:

    hat_P_n = zeros(N, 1)
    for i in 1:N
        hat_P_n[i] = prod((hat_P_nj[i, :]).^alpha_nj[i, :])
    end

    # Defining Lb:
    Lb = [sum(Ln[1:27]); Ln[28:69]]
    # Calculating the new equilibria of structures and land use:
    hat_H_n1=hat_H_n[1:27]
    hat_H_n2=zeros(27)
    try
        hat_H_n2 = eq4_h_endog(N,hat_H_n1,hat_omega_n0,phi_n,hat_P_n,hat_S_n,hat_L_n,beta_n,Lb)
        hat_H_n2[hat_H_n2 .< 0] = hat_H_n1[hat_H_n2 .< 0]
    catch
        hat_H_n2=hat_H_n1
    end
    hat_H_n1=hat_H_n2
    #hat_H_n1[ismissing.(hat_H_n1)].=1
    #hat_H_n1[hat_H_n1 .< 0].=1
    hat_H_n=[hat_H_n1;ones(N-27,1)]
    hat_H_n[isinf.(hat_H_n)].=1;hat_H_n[isnan.(hat_H_n)].=1;hat_H_n[iszero.(hat_H_n)].=1
    hat_H_n[hat_H_n .< 0.5].=1
    # Eq. 5:
    X_nk_l=eq5_H_var(N, J, alpha_nj, hat_omega_n0,hat_L_n,hat_H_n, beta_n, In,Ln, Sn, Sn_l,gamma_njk,pi_nij) #pi_nij_l

    # Eq 6: 
    #hat_omega_n=eq6_H_var(N, J, gamma_nj, pi_nij_l, X_nk_l, hat_L_n,hat_H_n, beta_n, In, Ln, Sn)
    num=zeros(N)
    den=zeros(N)
    hat_omega_n=zeros(N)
    for n in 1:N
        num[n] = sum(gamma_nj[n,j]*sum(pi_nij[i,n,j]*X_nk_l[i,j] for i=1:N) for j=1:J)
        den[n] = hat_H_n[n]^(beta_n[n])*hat_L_n[n]^(1 - beta_n[n]) * (In[n] * Ln[n] + Sn[n])
        hat_omega_n[n]=num[n]/den[n]
    end
    # 
    exc_dem = hat_omega_n0-hat_omega_n
    hat_omega_n1 = hat_omega_n0.*(1 .-0.05*exc_dem./hat_omega_n0).+(rand(Float64, (N))./100)
    epsilon=sum((hat_omega_n0-hat_omega_n1).^2)
    hat_omega_n0 = hat_omega_n1
    print(cont,"\n")
    print(epsilon,"\n")
    cont=cont+1
end

# Calculating TFP
tfp=zeros(N,J)
for n=1:N,j=1:J
    tfp[n,j]=exp(log((hat_T_nj[n,j]^gamma_nj[n,j])/((pi_nij_l[n,n,j]/pi_nij[n,n,j])^(1/theta_j[j]))))
end
hat_w_n=hat_omega_n0.*(hat_L_n.^(-beta_n))
DP=zeros(N*J,N)
PQ_vec=reshape(X_nk_l,1,J*N)'
for n=1:N
    DP[:,n]=reshape(pi_nij_l[n,:,:],N*J,1).*PQ_vec
end
aux_w=ones(J,1)*hat_w_n'
Exjn=zeros(J,N)
for j=1:J,n=1:N
    Exjn[j,n]=sum(DP[1+N*(j-1):N*J,n])'
end
PQ_vec0=reshape(X_nj',1,J*N)'
DPO=zeros(N*J,N)
for n=1:N
    DPO[:,n]=reshape(pi_nij[n,:,:],N*J,1).*PQ_vec0
end
Exjn0=zeros(J,N)
for n=1:N,j=1:J
    Exjn0[j,n]=sum(DPO[1+N*(j-1):N*j,n])
end
VALjn=gamma_nj.*(1 .- beta_n).*Exjn0'
hat_L_nj=((1 ./VALjn).*aux_w').*gamma_nj.*(1 .-beta_n).*Exjn'
# Calculating the GDP:
gdp=zeros(N,J)
for n=1:N,j=1:J
    gdp[n,j]=exp(log(tfp[n,j])+log(hat_L_nj[n,j])+log(hat_w_n[n]/hat_x_ij[n,j]))
end

# Aggregated Metrics
# Região
# TFP
tfp_n=zeros(N)
for n=1:N
    tfp_n[n]=sum(Y_nj[n,:]./sum(Y_nj[n,j] for j=1:J).*tfp[n,:])
end
# GDP
VARjn=(beta_n./(1 .-beta_n)).*VALjn
GDPnj_hat = ((hat_L_nj.*kron(ones(J,1),hat_omega_n0')')./hat_P_nj)
VAjn0  = VALjn + VARjn
VAj0 = sum(VAjn0',dims=2)
VAn0 = sum(VAjn0,dims=2)
VA0 = sum(VAj0)
GDP_n=zeros(N)
for n = 1:N
    GDP_n[n] = sum((VAjn0[n,:]./sum(VAjn0[n,j] for j=1:J)).*(GDPnj_hat[n,:]))
end


# Sector 
# TFP
tfp_j=zeros(J)
for j=1:J
    tfp_j[j]=sum(Y_nj[:,j]./sum(Y_nj[n,j] for n=1:N).*tfp[:,j])
end
tfp_j[isnan.(tfp_j)].=1
# GDP
GDPj=zeros(J)
for j = 1:J
    GDPj[j] = nansum((VAjn0[:,j]./sum(VAjn0[n,j] for n=1:N)).*(GDPnj_hat[:,j]))
end

# Aggregated GDP
AGDP = nansum((VAj0./VA0)'.*GDPj')
