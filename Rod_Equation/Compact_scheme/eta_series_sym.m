function [eta_ser_plus, eta_ser_minus] = eta_series_sym(w, nu, mu, N_eta)
%     nu = 20100; mu = 10000;
    a = 3/(12*mu + 4) - .5; b = 9*nu/(3*mu + 1) - 2; c = 1 - (12*nu + 3) / (6*mu + 2); d = 3*nu / (6*mu + 2);
%     syms w;

%     N_eta = 18;
    
    eta_leg = 0;
    eps_2 = (2*d-a*c)/a^2 + 2/a^2 * sqrt(-a*c*d + d^2 - 2*a^2*d^2 + a^2*b*d);
    w_1 = .5*(eps_2 - sqrt(eps_2^2-4));
    w_2 = .5*(eps_2 + sqrt(eps_2^2-4));
    for n = 0 : N_eta
        eta_leg = eta_leg + P(n, .5*(w_1 + w_2)) * w^n;
    end
    eta_ser_minus = 1/(2*d*w) * (-a*(1+w^2) - c*w - a*(1-w).*(w^2-(w_1+w_2)*w+1) * eta_leg); 
    eta_ser_plus = 1/(2*d*w) * (-a*(1+w^2) - c*w + a*(1-w).*(w^2-(w_1+w_2)*w+1) * eta_leg);
    eta_ser_plus = expand(eta_ser_plus);
    eta_ser_minus = expand(eta_ser_minus);
end
