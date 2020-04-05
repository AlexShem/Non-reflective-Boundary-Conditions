function [eta_ser_plus] = eta_series_num(nu, mu, N_eta)
%     a = 3/(12*mu + 4) - .5; b = 9*nu/(3*mu + 1) - 2; c = 1 - (12*nu + 3) / (6*mu + 2); d = 3*nu / (6*mu + 2);

    x = 1+72*nu/(6*mu-1)^2;
    leg_sum = zeros(1, N_eta + 2);
    for n = 0 : N_eta + 1
        leg_sum(end - n) = P(n, x);
    end
    eta_ser_plus = conv([1, -2*x, 1], leg_sum);
    eta_ser_plus = conv(eta_ser_plus, [-1, 1]);
    eta_ser_plus = polysum([1, -2*(-12*nu+6*mu-1)/(6*mu-1), 1], -eta_ser_plus);
    eta_ser_plus = (6*mu-1)/(12*nu) * eta_ser_plus;
    eta_ser_plus = eta_ser_plus(1 : end-1);
end
