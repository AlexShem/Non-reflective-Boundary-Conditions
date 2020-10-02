function [eta_ser_minus, eta_ser_plus] = series_eta(nu, mu, N_terms)
% eta_1 == eta minus (const + const * w + const * w^2 + ...)
% eta_2 == eta plus  (const + const * w + const * w^2 + ...) / w

x = 1+72*nu/(6*mu-1)^2;

legandre_sum = zeros(1, N_terms + 2);
for n = 0 : N_terms + 1
    legandre_sum(end - n) = P(n, x);
end

eta_ser = conv([1, -2*x, 1], legandre_sum);
eta_ser = conv(eta_ser, [-1, 1]);
eta_ser_minus = polysum([1, -2*(-12*nu+6*mu-1)/(6*mu-1), 1], -eta_ser);
eta_ser_plus = polysum([1, -2*(-12*nu+6*mu-1)/(6*mu-1), 1], eta_ser);
eta_ser_minus = (6*mu-1)/(12*nu) * eta_ser_minus;
eta_ser_plus = (6*mu-1)/(12*nu) * eta_ser_plus;

% Return polynomial of degree N_terms
% Note that the constant term inside of big brackets is zero
eta_ser_minus = eta_ser_minus(end - N_terms - 1 : end-1);
end
