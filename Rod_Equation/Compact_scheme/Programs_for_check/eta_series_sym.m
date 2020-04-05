function [eta_ser_plus, eta_ser_minus] = eta_series_sym(w, nu, mu, N_eta)
%     nu = 20100; mu = 10000;
%     a = 3/(12*mu + 4) - .5; b = 9*nu/(3*mu + 1) - 2; c = 1 - (12*nu + 3) / (6*mu + 2); d = 3*nu / (6*mu + 2);

    ep = 1 + 72 * nu / (6*mu - 1)^2;
    leg_sum = 0;
    for n = 0 : N_eta
        leg_sum = leg_sum + P(n, ep) * w^n;
    end
    eta_ser_minus = (6*mu - 1)/(12*nu*w) * ...
        (w^2 - 2*(-12*nu+6*mu-1)/(6*mu-1)*w + 1 ...
        - (1-w)*(w^2-2*ep*w+1) * leg_sum);
    eta_ser_plus = (6*mu - 1)/(12*nu*w) * ...
        (w^2 - 2*(-12*nu+6*mu-1)/(6*mu-1)*w + 1 ...
        + (1-w)*(w^2-2*ep*w+1) * leg_sum);
    eta_ser_minus = expand(eta_ser_minus);
    eta_ser_plus = expand(eta_ser_plus);
end
