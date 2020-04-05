function [lam_ser_1, lam_ser_2, lam_ser_3, lam_ser_4] = lambda_series_sym(w, nu, mu, N_lam, N_eta)
    a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*nu; s = nu/2;
    
    [eta_1, eta_2] = eta_series_sym(w, nu, mu, N_eta);
    eta_1_poly = sym2poly(eta_1);
    eta_2_poly = sym2poly(eta_2);
    
    r_13 = collect(w*poly2sym(eta_1_poly(1:end-1), w));
    r_24 = collect(w*poly2sym(eta_2_poly(1:end-1), w));
    
    % LAMBDA 1 3    -------------------------------------------------------
    theta_13 = 1/(2*s) * (-b - sqrt(b^2 - 4*s*(a-2*s)));
    mul = 1;
    mul_plus = 1;
    mul_minus = 1;
    sum_plus = 0;
    sum_minus = 0;

    for n = 0 : N_lam
        sum_plus = sum_plus + mul * mul_plus * r_13.^n;
        sum_minus = sum_minus + mul * mul_minus * r_13.^n;

        mul = -mul * (2*n+1) * (2*n+2) * (1-2*n) / ((-1 - 2*n) * (n+1)^2 * 4);
        mul_plus = mul_plus / (theta_13 + 2);
        mul_minus = mul_minus / (theta_13 - 2);
    end
    lam_ser_1 = eta_1/2 - sqrt(.25 * theta_13^2 - 1) * sum_plus .* sum_minus;
    lam_ser_3 = eta_1/2 + sqrt(.25 * theta_13^2 - 1) * sum_plus .* sum_minus;
    
    % LAMBDA 2 4    -------------------------------------------------------
    theta_24 = 1/(2*s) * (-b + sqrt(b^2 - 4*s*(a-2*s)));
    mul = 1;
    mul_plus = 1;
    mul_minus = 1;
    sum_plus = 0;
    sum_minus = 0;

    for n = 0 : N_lam
        sum_plus = sum_plus + mul * mul_plus * r_24.^n;
        sum_minus = sum_minus + mul * mul_minus * r_24.^n;

        mul = -mul * (2*n+1) * (2*n+2) * (1-2*n) / ((-1 - 2*n) * (n+1)^2 * 4);
        mul_plus = mul_plus / (theta_24 + 2);
        mul_minus = mul_minus / (theta_24 - 2);
    end
    lam_ser_2 = eta_2/2 - sqrt(.25 * theta_24^2 - 1) * sum_plus .* sum_minus;
    lam_ser_4 = eta_2/2 + sqrt(.25 * theta_24^2 - 1) * sum_plus .* sum_minus;
end