function [eta_ser_minus, eta_ser_plus] = eta_series_sym(w, nu, mu, N_eta)
    syms w;
    a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*nu; s = nu/2;

    eps_1 = 1/(b^2-4*s*(a-2*s)) * (-b*g+2*d*s + 2*sqrt(-b*g*d*s+d^2*s^2+a*g^2*s-2*g^2*s^2));
    eps_2 = 1/(b^2-4*s*(a-2*s)) * (-b*g+2*d*s - 2*sqrt(-b*g*d*s+d^2*s^2+a*g^2*s-2*g^2*s^2));
    
    geom_sum = 0;
    sum_1 = 0;
    sum_2 = 0;
    
    for n = 0 : N_eta
        geom_sum = geom_sum + (-1)^n * w.^(2*n);
        sum_1 = sum_1 + P(n, eps_1/2) * w.^n;
        sum_2 = sum_2 + P(n, eps_2/2) * w.^n;
    end

    eta_ser_minus = 1/(2*s) * geom_sum .* (-b*(1+w.^2) - g*w - sqrt(b^2 - 4*s*(a-2*s)) * (w.^2 - eps_1*w + 1) ...
        .* (w.^2 - eps_2*w + 1) .* sum_1 .* sum_2);
    eta_ser_plus = 1/(2*s) * geom_sum .* (-b*(1+w.^2) - g*w + sqrt(b^2 - 4*s*(a-2*s)) * (w.^2 - eps_1*w + 1) ...
        .* (w.^2 - eps_2*w + 1) .* sum_1 .* sum_2);
    
    eta_ser_minus = collect(eta_ser_minus);
    eta_ser_minus_poly = sym2poly(eta_ser_minus);
    eta_ser_minus = poly2sym(eta_ser_minus_poly(end-N_eta : end), w);
    
    eta_ser_plus = collect(eta_ser_plus);
    eta_ser_plus_poly = sym2poly(eta_ser_plus);
    eta_ser_plus = poly2sym(eta_ser_plus_poly(end-N_eta : end), w);
end