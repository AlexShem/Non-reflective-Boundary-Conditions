function [lam_1, lam_2, lam_3, lam_4] = lambda_series(w, nu, mu)
    a = 3/(12*mu + 4) - .5; b = 9*nu/(3*mu + 1) - 2; c = 1 - (12*nu + 3) / (6*mu + 2); d = 3*nu / (6*mu + 2);
    leg_deg = 30;
    [eta_1, eta_2] = eta_series(w, nu, mu, leg_deg);
    
%     Lambda 1, 3
    r = eta_1 + 1/a;
    lam_13_leg = 0;
    N_lam_13 = 20;
    for n = 0 : N_lam_13
        lam_13_leg = lam_13_leg + P(n, 1/sqrt(1-4*a^2)) * (a/sqrt(1-4*a^2))^n * r.^n;
    end
    lam_1 = eta_1/2 - sqrt(1/(4*a^2)-1) * (a^2/(1-4*a^2)*r.^2 - 2*a/(1-4*a^2)*r + 1) .* lam_13_leg;
    lam_3 = eta_1/2 + sqrt(1/(4*a^2)-1) * (a^2/(1-4*a^2)*r.^2 - 2*a/(1-4*a^2)*r + 1) .* lam_13_leg;
    
%     Lambda 2, 4
    N_frac = 20;
    N_lam_24 = 20;
    w_1 = (-b + 2*d - sqrt((b-2*d)^2-4))/2;
    w_2 = (-b + 2*d + sqrt((b-2*d)^2-4))/2;
    A = -1/sqrt((b-2*d)^2-4);
    fraction = 0;
    for k = 0 : N_frac
        fraction = fraction + (1/w_2^(k+1) - 1/w_1^(k+1)) * w.^k;
    end
    lam_2 = 0;
    mul_fac = 1;
    mul = 1;
    for n = 0 : N_lam_24
        lam_2 = lam_2 + mul * mul_fac / (1-2*n) * (w .* eta_1 .* fraction).^(2*n);
        mul = mul * (A*d)^2;
        mul_fac = mul_fac * (2*n+1)*(2*n+2)/(n+1)^2;
    end
    lam_4 = eta_2/2 .* (1 + lam_2);
    lam_2 = eta_2/2 .* (1 - lam_2);
end