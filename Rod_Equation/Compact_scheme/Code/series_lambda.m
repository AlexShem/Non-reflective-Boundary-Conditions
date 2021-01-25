function [lam_ser_1, lam_ser_2, lam_ser_3, lam_ser_4] = series_lambda(nu, mu, N_terms)
    a = 3/(12*mu + 4) - .5;
    b = (9*nu)/(3*mu + 1) - 2;
    c = 1 - (12*nu + 3)/(6*mu + 2);
    d = (3*nu)/(6*mu + 2);
    
    [eta_1, eta_2] = series_eta(nu, mu, N_terms + 5);
    
    r = [eta_1(1:end-1), 0]; % Remainder of \eta series
    
    % LAMBDA 1 3    -------------------------------------------------------
    mul = 1;
    mul_plus = 1;
    mul_minus = 1;
    sum_plus = zeros(1, N_terms+1); sum_plus(end) = mul * mul_plus;
    sum_minus = zeros(1, N_terms+1); sum_minus(end) = mul * mul_minus;

    r_13_cumul = r;
    for n = 0 : N_terms
        if n >= 2
            r_13_cumul = conv(r_13_cumul, r); r_13_cumul=r_13_cumul(end-N_terms:end);
            sum_plus = polysum(sum_plus , mul * mul_plus * r_13_cumul);
            sum_minus = polysum(sum_minus , mul * mul_minus * r_13_cumul);
        elseif n == 1
            sum_plus = polysum(sum_plus , mul * mul_plus * r);
            sum_minus = polysum(sum_minus , mul * mul_minus * r);
        end
        mul = mul * (2*n+1) * (2*n+2) * (1-2*n) / ((-1 - 2*n) * (n+1)^2 * 4);
        mul_plus = mul_plus / (1 + 2*a);
        mul_minus = mul_minus / (1 - 2*a);
    end
    lam_ser_1 = polysum(eta_1/2 , -1/(2*a) * sqrt(1 - 4*a^2) * conv(sum_plus , sum_minus));
    lam_ser_1 = lam_ser_1(end-N_terms:end);
    
    lam_ser_3 = polysum(eta_1/2 , +1/(2*a) * sqrt(1 - 4*a^2) * conv(sum_plus , sum_minus));
    lam_ser_3 = lam_ser_3(end-N_terms:end);
    
    % LAMBDA 2 4    -------------------------------------------------------
    mul = 1;
    mul_plus = 1;
    mul_minus = 1;
    sum_plus = zeros(1, N_terms+1); sum_plus(end) = mul * mul_plus;
    sum_minus = zeros(1, N_terms+1); sum_minus(end) = mul * mul_minus;
    
    A = (-3*mu - 1)/(2 * sqrt(3*nu * (3*nu - 6*mu - 2)));
    w1 = .5 * (-b + 2*d - sqrt((b - 2*d)^2 - 4));
    w2 = .5 * (-b + 2*d + sqrt((b - 2*d)^2 - 4));
    A = 1/(w1 - w2);

    eta_1_cumul = eta_1;
    inner_sum = flip(1:N_terms + 5); % +5 for excess (could be +0)
    inner_sum = w1.^inner_sum - w2.^inner_sum;
    inner_sum_conv = 1;
    
    for n = 0 : N_terms
        if n >= 2
            eta_1_cumul = conv(eta_1_cumul, eta_1); eta_1_cumul = eta_1_cumul(end-N_terms:end);
            inner_eta = conv(inner_sum, eta_1_cumul);
            sum_plus = polysum(sum_plus , mul_plus * inner_eta);
            sum_minus = polysum(sum_minus , mul_minus * inner_eta);
        elseif n == 1
            sum_plus = polysum(sum_plus , mul_plus * conv(inner_sum, eta_1));
            sum_minus = polysum(sum_minus , mul_minus * conv(inner_sum, eta_1));
        end
        mul = mul * A * d * (2*n+1) * (2*n+2) * (1-2*n) / ((-1 - 2*n) * (n+1)^2 * 2);
        mul_plus = mul;
        mul_minus = -mul;
        inner_sum_conv = conv(inner_sum_conv, inner_sum);
        % Reserve some memory and conv time
        inner_sum_conv = inner_sum_conv(end-length(inner_sum)+1 : end);
    end
    lam_ser_2 = .5 * conv(eta_2, polysum(1, -conv(sum_plus , sum_minus)));
    lam_ser_2 = lam_ser_2(end-N_terms:end);
    
    lam_ser_4 = .5 * conv(eta_2, polysum(1, +conv(sum_plus , sum_minus)));
    lam_ser_4 = lam_ser_4(end-N_terms:end);
end