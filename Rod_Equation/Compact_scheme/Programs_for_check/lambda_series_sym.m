function [lam_ser_1, lam_ser_2, lam_ser_3, lam_ser_4] = lambda_series_sym(w, nu, mu, N_lam, N_eta)
%     nu = 20100; mu = 10000;
%     nu = 40100; mu = 10000;
%     nu = 2.1/5; mu = 1/5;
%     nu = 9000; mu = 10000;
    a = 3/(12*mu + 4) - .5; b = 9*nu/(3*mu + 1) - 2; c = 1 - (12*nu + 3) / (6*mu + 2); d = 3*nu / (6*mu + 2);
%     syms w;
    
%     N_eta = 20;
    [eta_plus, eta_minus] = eta_series_sym(w, nu, mu, N_eta);
    
    eta_minus = sym2poly(eta_minus); %eta_minus = eta_minus(4:end); eta_minus = poly2sym(eta_minus, w);
%     eta_2 = sym2poly(eta_2); eta_2 = eta_2(4:end); eta_2 = poly2sym(eta_2, w);
%     r = eta_1 + 1/a;
%     eta_1_poly = sym2poly(eta_minus);
%     r = expand(w * poly2sym(eta_1_poly(1:end-1), w));
    eta_minus = eta_minus(end - N_eta + 1 : end);
    r = eta_minus(1 : end - 1);
%     N_lam_13_leg = 20;
    lam_13_leg = 0;
    for n = 0 : N_lam
        lam_13_leg = lam_13_leg + P(n, 1/sqrt(1-4*a^2)) * (a/sqrt(1-4*a^2))^n * r.^n;
    end
    lam_ser_1 = eta_minus/2 - sqrt(1/(4*a^2)-1) * (a^2/(1-4*a^2)*r.^2 - 2*a/(1-4*a^2)*r + 1) .* lam_13_leg;
    lam_ser_3 = eta_minus/2 + sqrt(1/(4*a^2)-1) * (a^2/(1-4*a^2)*r.^2 - 2*a/(1-4*a^2)*r + 1) .* lam_13_leg;
    
%     --------------------------------------------------------------------------------
    w_1 = (-b + 2*d - sqrt((b-2*d)^2-4))/2;
    w_2 = (-b + 2*d + sqrt((b-2*d)^2-4))/2;
    A = -1/sqrt((b-2*d)^2-4);
    fraction = 0;
%     N_frac = 10;
    for k = 0 : N_lam
        fraction = fraction + (1/w_2^(k+1) - 1/w_1^(k+1)) * w.^k;
    end
%     fraction = A * fraction;
    
%     N_lam_13_leg = 10;
    lam_24_leg = 0;
    mul_fac = 1;
    mul = 1;
    for n = 0 : N_lam
        lam_24_leg = lam_24_leg + mul * mul_fac / (1 - 2*n) * (w * eta_minus * fraction) ^ (2*n);
        mul = mul * (A*d)^2;
        mul_fac = mul_fac * (2*n + 1) * (2*n + 2) / (n + 1)^2;
    end
    lam_ser_2 = eta_plus / 2 * (1 - lam_24_leg);
    lam_ser_4 = eta_plus / 2 * (1 + lam_24_leg);
end