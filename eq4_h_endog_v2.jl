function eq4_h_endog(N,hat_H_n, hat_omega_n, phi_n, hat_P_n,  hat_S_n, hat_L_n, beta_n, Lb)
    num = zeros(27, 1)
    den = zeros(27, 1)
    hat_Ub=zeros(N-26)
    hat_Ub[1] = sum(Ln[1:27]./Lb[1].*(1 ./phi_n[1:27]).*(hat_omega_n[1:27].*hat_H_n[1:27].^(beta_n[1:27]).*hat_L_n[1:27].^(1 .-beta_n[1:27]))./hat_P_n[1:27].-((1 .-phi_n[1:27])./phi_n[1:27]).*hat_S_n[1:27]./hat_P_n[1:27])
    #hat_Ub[1]=1
    function equation_system!(F, hat_H_n)
    for n=1:27
        den[n] = (Lb[1]) * (hat_omega_n[n] / (phi_n[n] * hat_P_n[n] * hat_Ub[1] + (1 - phi_n[n]) * (hat_S_n[n] / hat_L_n[n]))) ^ (1 / beta_n[n])
        num[n] = hat_L_n[n] * sum(Ln[1:27] .*hat_H_n[1:27] .* ((hat_omega_n[1:27]) ./ (phi_n[1:27] .* hat_P_n[1:27] .* hat_Ub[1] .+ (1 .- phi_n[1:27]) .* (hat_S_n[1:27] ./ hat_L_n[1:27]))) .^ (1 ./ beta_n[1:27]))
        F[n]=hat_H_n[n]-num[n]/den[n]
    end
    end

    # Create an initial guess for the logarithm of the variables
    x0 = ones(27)

    # Solve the system of equations
    results = nlsolve(equation_system!, x0, ftol=1e-6)
    a=results.zero
    return a
end