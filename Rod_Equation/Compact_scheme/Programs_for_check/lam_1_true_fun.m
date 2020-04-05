function [lam_plus, lam_minus] = lam_1_true_fun(eta)
%     [eta_plus, eta_minus] = eta_true(nu, mu, omega)
    lam_plus = eta/2 + sqrt(eta.^2 / 4 - 1);
    lam_minus = eta/2 + sqrt(eta.^2 / 4 - 1);
end