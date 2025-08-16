# Main code
# Created by Alan Leal. E-mail: prof@alanleal-econ.com
# Setting random number generator seed
using Random
Random.seed!(12345)
# Setting main dimensions of analysis
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
# Scenario ii - Primary sectors surcharge, with trade cost reduction in all other sectors (do not run the two scenarios simultaneously)
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
using MAT


# Loading functions:
include("eq1.jl")
include("eq2.jl")
include("eq3.jl")
include("eq4_h_endog_v2.jl")
include("eq5_H_var.jl")
include("eq6_H_var.jl")
# Reading Julia files:
matfile = matopen(joinpath(@__DIR__, "land_emta_workspace.mat"))

# read all variables into a Dict
vars = read(matfile)  # this will give a Dict{String, Any}
close(matfile)

# assign each variable into Main (global environment)
for (name, value) in vars
    @eval Main $(Symbol(name)) = $value
end
theta_j=[theta_trad;0;0;0;0;0;0;0;0;0;0;0]
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
    Ln=[330718.29596203496; 269373.7540764827; 1.6342092793744425e6; 2.6409406239575744e6; 
    892485.6986513566; 204070.7939593827; 766550.8297753842; 1.2192035032174056e6; 6.2751926114648e6;
     3.982009057257481e6; 3.0868969110826906e6; 1.5575246580823185e6; 3.600423262728252e6; 
     1.3281988381358038e6; 959596.9739186086; 1.4053323082030565e6; 2.3776472906421227e6;
      2.9627233969011703e6; 2.2412426923678033e6; 1.3640796070238124e6; 1.611736216342898e6;
       8.92830198571447e6; 7.215148318221671e6; 2.307831066835844e7; 6.022835228605557e6;
        3.7742489898619703e6; 6.503607206113008e6; 9412.000600905001; 6092.612504000001; 
        6580.65336; 2961.945252; 16969.07724031; 4222.810013910001; 113085.44243957999; 
        585.496064; 312.92431923600003; 60792.57640000001; 582.474825; 31701.476120000003;
         965.6735467277122; 3803.63588; 40517.5984; 55702.83468; 7157.001048; 362.91449193999995;
          24.599652121600002; 12.017682402002999; 10949.566221823101; 3173.6775320000006; 
          37540.11088000001; 632.4439141727522; 20.09306803743; 2114.4589218529923; 
          518.1628400000001; 1587.3898323949802; 3443.0741831523; 247.66801199999998; 
          13222.492; 476.3852; 6677.0258176; 7514.434864000001; 3764.394924; 3236.2883481955;
           3213.0124076232005; 1479.7865880000002; 711.519828; 17184.278805392998; 
           609.7830377256719; 157084.392; 98368.59322589402]
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

# TFP
tfp=zeros(N,J)
for n=1:N,j=1:J
    tfp[n,j]=exp(log((hat_T_nj[n,j]^gamma_nj[n,j])/((pi_nij_l[n,n,j]/pi_nij[n,n,j])^(1/theta_j[j]))))
end
tfp_n=zeros(N)
for n=1:N
    tfp_n[n]=sum(Y_nj[n,:]./sum(Y_nj[n,j] for j=1:J).*tfp[n,:])
end
# Relevant Metrics
welfare=zeros(N)
for n=1:N
    wn=hat_omega_n0[n]
    welfare[n]=exp(sum(alpha_nj[n,:].*(log.(tfp[n,:])-log.(wn./hat_x_ij[n,:]))))
end
hat_H_n
