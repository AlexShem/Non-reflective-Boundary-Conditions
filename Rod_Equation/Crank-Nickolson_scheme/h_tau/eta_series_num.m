function [eta_ser_minus, eta_ser_plus] = eta_series_num(nu, mu, N_eta)
%     a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*mu; s = nu/2;
%     eps_1 = 1/(b^2-4*s*(a-2*s)) * (-b*g+2*d*s + 2*sqrt(-b*g*d*s+d^2*s^2+a*g^2*s-2*g^2*s^2));
    eps_1 = mu.^2./(mu.^2 - 2*nu);

    geom_sum = zeros(1, 2*N_eta + 1);
    sum_1 = zeros(1, N_eta+1);
    
    for n = 0 : N_eta
        geom_sum(end - 2*n) = (-1)^n;
        sum_1(end-n) = P(n, eps_1);
    end
    
    % —юда можно приделать обрезание радов
    conv_sums = conv([1, -2*eps_1, 1], sum_1);
    if length(conv_sums) > N_eta+1
        conv_sums = conv_sums(end-N_eta:end);
    end
%     eta_ser_minus = polysum([-b,-g,-b], -sqrt(b^2 - 4*s*(a-2*s)) * conv_sums);
    eta_ser_minus = polysum(conv(conv_sums, -sqrt(mu^2-2*nu)*[-1 1]), [mu+2*nu, -2*mu, mu+2*nu]);
    eta_ser_minus = 1/nu * conv(geom_sum,  eta_ser_minus);
    if length(eta_ser_minus) > N_eta+1
        eta_ser_minus = eta_ser_minus(end-N_eta:end);
    end
    
    eta_ser_plus = polysum(conv(conv_sums, sqrt(mu^2-2*nu)*[-1 1]), [mu+2*nu, -2*mu, mu+2*nu]);
    eta_ser_plus = 1/nu * conv(geom_sum,  eta_ser_plus);
    if length(eta_ser_plus) > N_eta+1
        eta_ser_plus = eta_ser_plus(end-N_eta:end);
    end
end