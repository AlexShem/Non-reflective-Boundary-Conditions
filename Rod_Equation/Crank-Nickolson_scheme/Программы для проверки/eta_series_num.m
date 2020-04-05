function [eta_ser_minus, eta_ser_plus] = eta_series_num(nu, mu, N_eta)
    a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*nu; s = nu/2;
    eps_1 = 1/(b^2-4*s*(a-2*s)) * (-b*g+2*d*s + 2*sqrt(-b*g*d*s+d^2*s^2+a*g^2*s-2*g^2*s^2));
    eps_2 = 1/(b^2-4*s*(a-2*s)) * (-b*g+2*d*s - 2*sqrt(-b*g*d*s+d^2*s^2+a*g^2*s-2*g^2*s^2));
    
    geom_sum = zeros(1, 2*N_eta + 1);
    sum_1 = zeros(1, N_eta+1);
    sum_2 = zeros(1, N_eta+1);
    
    for n = 0 : N_eta
        geom_sum(end - 2*n) = (-1)^n;
        sum_1(end-n) = P(n, eps_1/2);
        sum_2(end-n) = P(n, eps_2/2);
    end
    
    % —юда можно приделать обрезание р€дов
    conv_sums = conv(conv([1, -eps_1, 1], [1, -eps_2, 1]), conv(sum_1, sum_2));
    if length(conv_sums) > N_eta+1
        conv_sums = conv_sums(end-N_eta:end);
    end
    eta_ser_minus = polysum([-b,-g,-b], -sqrt(b^2 - 4*s*(a-2*s)) * conv_sums);
    eta_ser_minus = 1/(2*s) * conv(geom_sum,  eta_ser_minus);
    if length(eta_ser_minus) > N_eta+1
        eta_ser_minus = eta_ser_minus(end-N_eta:end);
    end

    eta_ser_plus = polysum([-b,-g,-b], sqrt(b^2 - 4*s*(a-2*s)) * conv_sums);
    eta_ser_plus = 1/(2*s) * conv(geom_sum,  eta_ser_plus);
    if length(eta_ser_plus) > N_eta+1
        eta_ser_plus = eta_ser_plus(end-N_eta:end);
    end
end