function [eta_1, eta_2] = eta_series(w, nu, mu, leg_deg)
    a = 3/(12*mu + 4) - .5; b = 9*nu/(3*mu + 1) - 2; c = 1 - (12*nu + 3) / (6*mu + 2); d = 3*nu / (6*mu + 2);
    eps_1 = (-2*a*c + 4*d)/a^2 - 2;
    w_1 = .5*(eps_1 - sqrt(eps_1^2-4));
    w_2 = .5*(eps_1 + sqrt(eps_1^2-4));
    eta_2 = zeros(size(w));
    for n = 0 : leg_deg
        eta_2 = eta_2 + P(n, .5*eps_1) * w.^n;
    end
    eta_1 = 1./(2*d*w) .* (-a*(1+w.^2) - c*w + a*(1-w).*(w.^2-eps_1*w+1).*eta_2);
    eta_2 = 1./(2*d*w) .* (-a*(1+w.^2) - c*w - a*(1-w).*(w.^2-eps_1*w+1).*eta_2); 
end