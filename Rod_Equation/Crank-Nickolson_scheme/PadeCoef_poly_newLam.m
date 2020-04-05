function [polyCoeff_bor, polyCoeff_prebor] = PadeCoef_poly_newLam(nu, mu, deg_bor, deg_prebor)
    step_num_1 = length(deg_bor);
    step_num_2 = length(deg_prebor);
    totalDeg_1 = sum(deg_bor + 1);
    totalDeg_2 = sum(deg_prebor + 1);
    maxTotalDeg = max([totalDeg_1, totalDeg_2]);
    
    [~, ~, lam_ser_3, lam_ser_4, sqrt_lam13, sqrt_lam24] = lambda_series_num_newLam(nu, mu, maxTotalDeg+5);
    [eta_ser_1, eta_ser_2] = eta_series_num(nu, mu, maxTotalDeg+5);
    
    num_unknCoeff_bor = totalDeg_1;
    num_unknCoeff_prebor = totalDeg_2;
    
    A_bor = zeros(num_unknCoeff_bor);
    for i = 0 : step_num_1 - 1
        lam_3_conv = lam_ser_3;
        lam_4_conv = lam_ser_4;
        if i == 2
            lam_3_conv = polysum(-1, polysum(conv(eta_ser_1, eta_ser_1)/2, conv(eta_ser_1, sqrt_lam13)));
            lam_4_conv = polysum(-1, polysum(conv(eta_ser_2, eta_ser_2)/2, conv(eta_ser_2, sqrt_lam24)));
        elseif i == 3
            lam_3_conv = polysum(polysum(conv(conv(eta_ser_1, eta_ser_1), eta_ser_1)/2, -eta_ser_1*1.5), ...
                conv(polysum(conv(eta_ser_1, eta_ser_1), -1), sqrt_lam13));
            lam_4_conv = polysum(polysum(conv(conv(eta_ser_2, eta_ser_2), eta_ser_2)/2, -eta_ser_2*1.5), ...
                conv(polysum(conv(eta_ser_2, eta_ser_2), -1), sqrt_lam24));
        else
            for j = 2 : i
                lam_3_conv = conv(lam_3_conv, lam_ser_3);
                lam_4_conv = conv(lam_4_conv, lam_ser_4);
            end
        end
        for k = 0 : deg_bor(i+1)
            subs_var = zeros(1, deg_bor(i+1)+1);
            subs_var(end-k) = 1;
            if i ~= 0
                tmp_conv_3 = conv(lam_3_conv, subs_var);
                tmp_conv_4 = conv(lam_4_conv, subs_var);
            else
                tmp_conv_3 = subs_var;
                tmp_conv_4 = subs_var;
            end
            if length(tmp_conv_3) > num_unknCoeff_bor
                tmp_conv_3 = tmp_conv_3(end - num_unknCoeff_bor + 1 : end);
            elseif length(tmp_conv_3) < num_unknCoeff_bor
                tmp_conv_3 = [zeros(1, num_unknCoeff_bor - length(tmp_conv_3)), tmp_conv_3];
            end
            if length(tmp_conv_4) > num_unknCoeff_bor
                tmp_conv_4 = tmp_conv_4(end - num_unknCoeff_bor + 1 : end);
%                 tmp_conv_4 = tmp_conv_4(end - num_unknCoeff_bor-1 : end-2);
            elseif length(tmp_conv_4) < num_unknCoeff_bor
                tmp_conv_4 = [zeros(1, num_unknCoeff_bor - length(tmp_conv_4)), tmp_conv_4];
            end
%             A_bor(:, sum(deg_bor(1:i)+1) + k+1) = flip(tmp_conv);
            one_bound = ceil(num_unknCoeff_bor/2 - 1);
            two_bound = floor(num_unknCoeff_bor/2 - 1);
            A_bor(3:3+one_bound-1, sum(deg_bor(1:i)+1) + k+1) = ...
                (tmp_conv_3(end:-1:end-one_bound+1));
            A_bor(3+one_bound:end, sum(deg_bor(1:i)+1) + k+1) = ...
                (tmp_conv_4(end:-1:end-two_bound+1));
        end
    end
%     A_bor(1:2, :) = 0;
    A_bor(1, 1) = 1; A_bor(2, deg_bor(1)+2) = 1;
    tmp = A_bor;
    A_bor(3:3+one_bound-1, :) = (tmp(3:3+one_bound-1, :) + tmp(3+one_bound:end, :)) / 2;
    A_bor(3+one_bound:end, :) = (tmp(3:3+one_bound-1, :) - tmp(3+one_bound:end, :)) / 2i;
    clear tmp;
    B_bor = zeros(num_unknCoeff_bor, 1);
    B_bor(1, 1) = 1;
    
    A_prebor = zeros(num_unknCoeff_prebor);
    for i = 0 : step_num_2 - 1
        lam_3_conv = lam_ser_3;
        lam_4_conv = lam_ser_4;
        for j = 2 : i
            lam_3_conv = conv(lam_3_conv, lam_ser_3);
            lam_4_conv = conv(lam_4_conv, lam_ser_4);
        end
        for k = 0 : deg_prebor(i+1)
            subs_var = zeros(1, deg_prebor(i+1)+1);
            subs_var(end-k) = 1;
            if i ~= 0
                tmp_conv_3 = conv(lam_3_conv, subs_var);
                tmp_conv_4 = conv(lam_4_conv, subs_var);
            else
                tmp_conv_3 = subs_var;
                tmp_conv_4 = subs_var;
            end
            if length(tmp_conv_3) > num_unknCoeff_prebor
                tmp_conv_3 = tmp_conv_3(end - num_unknCoeff_prebor + 1 : end);
            elseif length(tmp_conv_3) < num_unknCoeff_prebor
                tmp_conv_3 = [zeros(1, num_unknCoeff_prebor - length(tmp_conv_3)), tmp_conv_3];
            end
            if length(tmp_conv_4) > num_unknCoeff_prebor
                tmp_conv_4 = tmp_conv_4(end - num_unknCoeff_prebor + 1 : end);
%                 tmp_conv_4 = tmp_conv_4(end - num_unknCoeff_bor-1 : end-2);
            elseif length(tmp_conv_4) < num_unknCoeff_prebor
                tmp_conv_4 = [zeros(1, num_unknCoeff_prebor - length(tmp_conv_4)), tmp_conv_4];
            end
%             A_bor(:, sum(deg_bor(1:i)+1) + k+1) = flip(tmp_conv);
            one_bound = ceil(num_unknCoeff_prebor/2 - 1);
            two_bound = floor(num_unknCoeff_prebor/2 - 1);
            A_prebor(3:3+one_bound-1, sum(deg_prebor(1:i)+1) + k+1) = ...
                (tmp_conv_3(end:-1:end-one_bound+1));
            A_prebor(3+one_bound:end, sum(deg_prebor(1:i)+1) + k+1) = ...
                (tmp_conv_4(end:-1:end-two_bound+1));
        end
    end
%     A_prebor(1:2, :) = 0;
    A_prebor(1, deg_prebor(1)+2) = 1; A_prebor(2, 1) = 1;
    tmp = A_prebor;
    A_prebor(3:3+one_bound-1, :) = (tmp(3:3+one_bound-1, :) + tmp(3+one_bound:end, :)) / 2;
    A_prebor(3+one_bound:end, :) = (tmp(3:3+one_bound-1, :) - tmp(3+one_bound:end, :)) / 2i;
    clear tmp;
    B_prebor = zeros(num_unknCoeff_prebor, 1);
    B_prebor(1, 1) = 1;
    
    % A bit of cheating
    A_bor = real(A_bor); A_prebor = real(A_prebor);
    
    eig_val_bor = eig(A_bor); eig_bor_min = min(abs(eig_val_bor));
    eig_val_prebor = eig(A_prebor); eig_prebor_min = min(abs(eig_val_prebor));

    if eig_bor_min <= 1e-10 && eig_prebor_min <= 1e-10
        error('Both sets of border points have no Pade approximation');
    elseif eig_bor_min <= 1e-10
        error('Border point set has no Pade aproximation');
    elseif eig_prebor_min <= 1e-10
        error('Preborder point set has no Pade aproximation');
    end
    
    polyCoeff_bor = linsolve(A_bor, B_bor);
    polyCoeff_prebor = linsolve(A_prebor, B_prebor);
end