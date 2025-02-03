# Coefs

using MAT
using LinearAlgebra
using Statistics
# Changing the working directory:
cd("/Users/alanleal/Library/Mobile Documents/com~apple~CloudDocs/Tese/Modelos e Códigos")
# Reading the data:
data_mat=matread("matrix_datav2.mat")
io_mat=data_mat["io"]
sa_mat=data_mat["sa"]
# Defining some objects of interest:
N=70
J=26

# Creating the objects of interest:
# 2.1) X_nij
# The dimension of this array is: 70 regions by 70 regions by 26 sectors
X_nij = zeros(N, N, 26)
for x in 1:N
    for z in 1:N
        a = io_mat[(26*x-25):(26*x), (26*z-25):(26*z)]
        for w in 1:26
            X_nij[z, x, w] = sum(a[:, w])
        end
    end
end

# 2.2) pi_nij
# Now we can calculate the share of region n's expenditure on sector j from region i.
pi_nij = zeros(N, N, 26)
for w in 1:26
    for z in 1:N
        a = X_nij[z, :, w]
        for x in 1:N
            pi_nij[z, x, w] = X_nij[z, x, w] / sum(a)
        end
    end
end
pi_nij[isnan.(pi_nij)].=0


# 2.3) S_n (trade balance)
Sn = zeros(N, 1)
for i in 1:N
    a = X_nij[i, :, :]
    b = X_nij[:, i, :]
    Sn[i] = sum(a) - sum(b)
end

# 2.4) In (per capita income)
# Per capita income is given by: In = (VAn - Sn) / Ln
In = zeros(N, 1)
Vn = zeros(N, 1)
Ln = zeros(N, 1)
for i in 1:N
    Vn[i] = sum(sa_mat[4, (26*i-25):(26*i)])
    Ln[i] = sum(sa_mat[1, (26*i-25):(26*i)])
    In[i] = (Vn[i] - Sn[i]) / Ln[i]
end

# Mudança PROVISÓRIA no código: pesquisar como consertar problema para Russia, Slovakia e Slovenia
in_media=mean(In[28:70])
In[In.<0].=in_media
# PROVISÖRIO: PENSAR EM COMO CONSERTAR ESSA RENDA PER CAPITA NEGATIVA

Y_nj=zeros(N,J)
for n in 1:N,j in 1:J
    Y_nj[n,j]=sa_mat[4,26*(n-1)+j]
end

# 2.5) beta_n
beta_n = zeros(N, 1)
wl = zeros(N, 1)
for i in 1:N
    wl[i] = sum(sa_mat[3, (26*i-25):(26*i)])
    beta_n[i] = wl[i] / Vn[i]
end

# 2.6) gamma_nj
VA_nj = zeros(N, J)
PIB_nj = zeros(N, J)
gamma_nj = zeros(N, J)
for i in 1:N
    VA_nj[i, :] = sa_mat[4, (26*i-25):(26*i)]
    PIB_nj[i, :] = sa_mat[2, (26*i-25):(26*i)]
    gamma_nj[i, :] = VA_nj[i, :] ./ PIB_nj[i, :]
end
gamma_nj[isnan.(gamma_nj)].=0

# 2.7) gamma_njk
# Creating the M matrix first:
M_njk = zeros(N, J, J)
for n in 1:N
    for j in 1:J
        for k in 1:J
            M_njk[n, j, k] = io_mat[(26*n-26+k), (26*n-26+j)]
        end
    end
end


gamma_njk = zeros(N, J, J)
for z in 1:26
    for l in 1:N
        for w in 1:26
            gamma_njk[l, w, z] = M_njk[l, w, z] / PIB_nj[l, w]
        end
    end
end
gamma_njk[isnan.(gamma_njk)].=0


# 2.8) X_nj
X_nj = zeros(N, J)
for i in 1:N,j in 1:J
        X_nj[i, j] = sum(X_nij[i,:,j])
end


# 2.9) alpha_nj:
alpha_nj = zeros(N, J)
for n in 1:N,j in 1:J
        alpha_nj[n, j] = (X_nj[n,j]-sum(gamma_njk[n,:,j].*(sum(pi_nij[:,n,:].*X_nj,dims=1)')))/(In[n]*Ln[n])
end

alpha_nj[isnan.(alpha_nj)].=0
alpha_nj[isinf.(alpha_nj)].=0
alpha_nj[alpha_nj.<0].=0

# 2.10) theta_j (dispersão de produtividade entre as firmas)
# Vou considerar as estimativas do Maggi, enquanto eu não as estimo eu mesmo
# Já estimei, preciso validar com o Haddad
theta_trad=[6.092008;14.222992;5.409076;5.430116;5.375593;7.03864;15.829086;6.901034;4.295246;11.370322;5.907845;42.535307;3.091325;9.291044;15.667740]
theta_j = [theta_trad;ones(J-15,1)]
theta_j=vec(theta_j)
# Agora, vamos criar o vetor da dispersão de produtividade corretamente. 

# 2.11) phi_n
phi_n = zeros(N, 1)
for n in 1:N
    phi_n[n] = 1 / (1 + Sn[n] / (In[n] * Ln[n]))
end
