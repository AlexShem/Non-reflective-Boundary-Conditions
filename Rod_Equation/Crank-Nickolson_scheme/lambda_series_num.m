function [lam_ser_1, lam_ser_2, lam_ser_3, lam_ser_4] = lambda_series_num(nu, mu, N_lam)
    a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*mu; s = nu/2;
    
    [eta_1, eta_2] = eta_series_num(nu, mu, N_lam);
    
    r_13 = [eta_1(1:end-1), 0];
    r_24 = [eta_2(1:end-1), 0];
    
    % LAMBDA 1 3    -------------------------------------------------------
    theta_13 = 1/(2*s) * (-b - sqrt(b^2 - 4*s*(a-2*s)));
    mul = 1;
    mul_plus = 1;
    mul_minus = 1;
    sum_plus = zeros(1, N_lam+1); sum_plus(end) = mul * mul_plus;
    sum_minus = zeros(1, N_lam+1); sum_minus(end) = mul * mul_minus;

    r13_cumul = r_13;
    for n = 0 : N_lam
        if n >= 2
            r13_cumul = conv(r13_cumul, r_13); r13_cumul=r13_cumul(end-N_lam:end);
            sum_plus = polysum(sum_plus , mul * mul_plus * r13_cumul);
            sum_minus = polysum(sum_minus , mul * mul_minus * r13_cumul);
        elseif n == 1
            sum_plus = polysum(sum_plus , mul * mul_plus * r_13);
            sum_minus = polysum(sum_minus , mul * mul_minus * r_13);
        end
        mul = -mul * (2*n+1) * (2*n+2) * (1-2*n) / ((-1 - 2*n) * (n+1)^2 * 4);
        mul_plus = mul_plus / (theta_13 + 2);
        mul_minus = mul_minus / (theta_13 - 2);
    end
    lam_ser_1 = polysum(eta_1/2 , -sqrt(.25 * theta_13^2 - 1) * conv(sum_plus , sum_minus));
    lam_ser_1 = lam_ser_1(end-N_lam:end);
    lam_ser_3 = polysum(eta_1/2 , sqrt(.25 * theta_13^2 - 1) * conv(sum_plus , sum_minus));
    lam_ser_3 = lam_ser_3(end-N_lam:end);
    
    % LAMBDA 2 4    -------------------------------------------------------
    theta_24 = 1/(2*s) * (-b + sqrt(b^2 - 4*s*(a-2*s)));
    mul = 1;
    mul_plus = 1;
    mul_minus = 1;
    sum_plus = zeros(1, N_lam+1); sum_plus(end) = mul * mul_plus;
    sum_minus = zeros(1, N_lam+1); sum_minus(end) = mul * mul_minus;

    r_24_cumul = r_24;
    for n = 0 : N_lam
        if n >= 2
            r_24_cumul = conv(r_24_cumul, r_24); r_24_cumul = r_24_cumul(end-N_lam:end);
            sum_plus = polysum(sum_plus , mul * mul_plus * r_24_cumul);
            sum_minus = polysum(sum_minus , mul * mul_minus * r_24_cumul);
        elseif n == 1
            sum_plus = polysum(sum_plus , mul * mul_plus * r_24);
            sum_minus = polysum(sum_minus , mul * mul_minus * r_24);
        end
        mul = -mul * (2*n+1) * (2*n+2) * (1-2*n) / ((-1 - 2*n) * (n+1)^2 * 4);
        mul_plus = mul_plus / (theta_24 + 2);
        mul_minus = mul_minus / (theta_24 - 2);
    end
    lam_ser_2 = polysum(eta_2/2 , -sqrt(.25 * theta_24^2 - 1) * conv(sum_plus , sum_minus));
    lam_ser_2 = lam_ser_2(end-N_lam:end);
    lam_ser_4 = polysum(eta_2/2 , sqrt(.25 * theta_24^2 - 1) * conv(sum_plus , sum_minus));
    lam_ser_4 = lam_ser_4(end-N_lam:end);
end