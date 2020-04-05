function [polyCoeff_bor, polyCoeff_prebor] = ...
    PadeCoef_poly_const_x_1w2(nu, mu, deg_bor, deg_prebor, constant_condition, x_condition)

    step_num_1 = length(deg_bor);
    step_num_2 = length(deg_prebor);
    totalDeg_1 = sum(deg_bor + 1);
    totalDeg_2 = sum(deg_prebor + 1);
    maxTotalDeg = max([totalDeg_1, totalDeg_2]);
    
    % Use condition for u(t, x) = const
	condition_type = xor(constant_condition, x_condition);
    if condition_type == 1 && (mod(totalDeg_1, 2) == 0)
        error('Number of coefs in Border Deg must be odd');
    elseif condition_type == 1 && (mod(totalDeg_2, 2) == 0)
        error('Number of coefs in PreBorder Deg must be odd');
    elseif condition_type == 0 && (mod(totalDeg_1, 2) == 1)
        error('Number of coefs in Border Deg must be even');
    elseif condition_type == 0 && (mod(totalDeg_2, 2) == 1)
        error('Number of coefs in PreBorder Deg must be even');
    end
    
    [~, ~, lam_ser_3, lam_ser_4] = lambda_series_num(nu, mu, maxTotalDeg+5);

    num_unknCoeff_bor = totalDeg_1;
    num_unknCoeff_prebor = totalDeg_2;
    
    A_bor = zeros(num_unknCoeff_bor);
    for i = 0 : step_num_1 - 1
        lam_3_conv = conv(lam_ser_3, [1 0 1]);
%         lam_4_conv = conv(lam_ser_4, [1 0 1]);
        lam_4_conv = conv(lam_ser_4, [1 0 1]);
        for j = 2 : i
            lam_3_conv = conv(conv(lam_3_conv, lam_ser_3), [1 0 1]);
            lam_4_conv = conv(conv(lam_4_conv, lam_ser_4), [1 0 1]);
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
            elseif length(tmp_conv_4) < num_unknCoeff_bor
                tmp_conv_4 = [zeros(1, num_unknCoeff_bor - length(tmp_conv_4)), tmp_conv_4];
            end
            
            switch condition_type
                case 0
                    if x_condition == 0
                        one_bound = num_unknCoeff_bor/2 - 1;
                        two_bound = num_unknCoeff_bor/2 - 1;
                        A_bor(3:3+one_bound-1, sum(deg_bor(1:i)+1) + k+1) = ...
                            (tmp_conv_3(end:-1:end-one_bound+1));
                        A_bor(3+one_bound:end, sum(deg_bor(1:i)+1) + k+1) = ...
                            (tmp_conv_4(end:-1:end-two_bound+1));
                    else
                        one_bound = num_unknCoeff_bor/2 - 2;
                        two_bound = num_unknCoeff_bor/2 - 2;
                        A_bor(5:5+one_bound-1, sum(deg_bor(1:i)+1) + k+1) = ...
                            (tmp_conv_3(end:-1:end-one_bound+1));
                        A_bor(5+one_bound:end, sum(deg_bor(1:i)+1) + k+1) = ...
                            (tmp_conv_4(end:-1:end-two_bound+1));
                    end
                case 1
                    one_bound = (num_unknCoeff_bor - 1)/2 - 1;
                    two_bound = (num_unknCoeff_bor - 1)/2 - 1;
                    A_bor(4:4+one_bound-1, sum(deg_bor(1:i)+1) + k+1) = ...
                        (tmp_conv_3(end:-1:end-one_bound+1));
                    A_bor(4+one_bound:end, sum(deg_bor(1:i)+1) + k+1) = ...
                        (tmp_conv_4(end:-1:end-two_bound+1));
            end
        end
    end
    % Normalization of coefficients
    A_bor(1, 1) = 1; A_bor(2, deg_bor(1)+2) = 1;
    % Constant condition u(t, x) = const
    tmp = A_bor;
    if constant_condition == 0 && x_condition == 0
        A_bor(3:3+one_bound-1, :) = (tmp(3:3+one_bound-1, :) + tmp(3+one_bound:end, :)) / 2;
        A_bor(3+one_bound:end, :) = (tmp(3:3+one_bound-1, :) - tmp(3+one_bound:end, :)) / 2i;
    elseif constant_condition == 1 && x_condition == 0
        A_bor(4:4+one_bound-1, :) = (tmp(4:4+one_bound-1, :) + tmp(4+one_bound:end, :)) / 2;
        A_bor(4+one_bound:end, :) = (tmp(4:4+one_bound-1, :) - tmp(4+one_bound:end, :)) / 2i;
        A_bor(3,:) = 1;
    elseif constant_condition == 0 && x_condition == 1
        A_bor(4:4+one_bound-1, :) = (tmp(4:4+one_bound-1, :) + tmp(4+one_bound:end, :)) / 2;
        A_bor(4+one_bound:end, :) = (tmp(4:4+one_bound-1, :) - tmp(4+one_bound:end, :)) / 2i;
        for k = 1 : step_num_1
            A_bor(3, sum(deg_bor(1:(k-1))+1) + 1 : sum(deg_bor(1:k)+1)) = k-1;
        end
    elseif constant_condition == 1 && x_condition == 1
        A_bor(5:5+one_bound-1, :) = (tmp(5:5+one_bound-1, :) + tmp(5+one_bound:end, :)) / 2;
        A_bor(5+one_bound:end, :) = (tmp(5:5+one_bound-1, :) - tmp(5+one_bound:end, :)) / 2i;
        for k = 1 : step_num_1
            A_bor(4, sum(deg_bor(1:(k-1))+1) + 1 : sum(deg_bor(1:k)+1)) = k-1;
        end
        A_bor(3, :) = 1;
    end
    clear tmp;
    B_bor = zeros(num_unknCoeff_bor, 1);
    B_bor(1, 1) = 1;
    
% %     ----------------------------------------------------------------------------------------------
    A_prebor = zeros(num_unknCoeff_prebor);
    for i = 0 : step_num_2 - 1
        lam_3_conv = conv(lam_ser_3, [1 0 1]);
        lam_4_conv = conv(lam_ser_4, [1 0 1]);
        for j = 2 : i
            lam_3_conv = conv(conv(lam_3_conv, lam_ser_3), [1 0 1]);
            lam_4_conv = conv(conv(lam_4_conv, lam_ser_4), [1 0 1]);
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
            elseif length(tmp_conv_4) < num_unknCoeff_prebor
                tmp_conv_4 = [zeros(1, num_unknCoeff_prebor - length(tmp_conv_4)), tmp_conv_4];
            end
            
            switch condition_type
                case 0
                    if x_condition == 0
                        one_bound = num_unknCoeff_bor/2 - 1;
                        two_bound = num_unknCoeff_bor/2 - 1;
                        A_prebor(3:3+one_bound-1, sum(deg_bor(1:i)+1) + k+1) = ...
                            (tmp_conv_3(end:-1:end-one_bound+1));
                        A_prebor(3+one_bound:end, sum(deg_bor(1:i)+1) + k+1) = ...
                            (tmp_conv_4(end:-1:end-two_bound+1));
                    else
                        one_bound = num_unknCoeff_bor/2 - 2;
                        two_bound = num_unknCoeff_bor/2 - 2;
                        A_prebor(5:5+one_bound-1, sum(deg_bor(1:i)+1) + k+1) = ...
                            (tmp_conv_3(end:-1:end-one_bound+1));
                        A_prebor(5+one_bound:end, sum(deg_bor(1:i)+1) + k+1) = ...
                            (tmp_conv_4(end:-1:end-two_bound+1));
                    end
                case 1
                    one_bound = (num_unknCoeff_bor - 1)/2 - 1;
                    two_bound = (num_unknCoeff_bor - 1)/2 - 1;
                    A_prebor(4:4+one_bound-1, sum(deg_bor(1:i)+1) + k+1) = ...
                        (tmp_conv_3(end:-1:end-one_bound+1));
                    A_prebor(4+one_bound:end, sum(deg_bor(1:i)+1) + k+1) = ...
                        (tmp_conv_4(end:-1:end-two_bound+1));
            end
        end
    end
%     A_prebor(1:2, :) = 0;
    A_prebor(1, deg_prebor(1)+2) = 1; A_prebor(2, 1) = 1;
    tmp = A_prebor;
    if constant_condition == 0 && x_condition == 0
        A_prebor(3:3+one_bound-1, :) = (tmp(3:3+one_bound-1, :) + tmp(3+one_bound:end, :)) / 2;
        A_prebor(3+one_bound:end, :) = (tmp(3:3+one_bound-1, :) - tmp(3+one_bound:end, :)) / 2i;
    elseif constant_condition == 1 && x_condition == 0
        A_prebor(4:4+one_bound-1, :) = (tmp(4:4+one_bound-1, :) + tmp(4+one_bound:end, :)) / 2;
        A_prebor(4+one_bound:end, :) = (tmp(4:4+one_bound-1, :) - tmp(4+one_bound:end, :)) / 2i;
        A_prebor(3,:) = 1;
    elseif constant_condition == 0 && x_condition == 1
        A_prebor(4:4+one_bound-1, :) = (tmp(4:4+one_bound-1, :) + tmp(4+one_bound:end, :)) / 2;
        A_prebor(4+one_bound:end, :) = (tmp(4:4+one_bound-1, :) - tmp(4+one_bound:end, :)) / 2i;
        for k = 1 : step_num_1
            A_prebor(3, sum(deg_bor(1:(k-1))+1) + 1 : sum(deg_bor(1:k)+1)) = k-1;
        end
    elseif constant_condition == 1 && x_condition == 1
        A_prebor(5:5+one_bound-1, :) = (tmp(5:5+one_bound-1, :) + tmp(5+one_bound:end, :)) / 2;
        A_prebor(5+one_bound:end, :) = (tmp(5:5+one_bound-1, :) - tmp(5+one_bound:end, :)) / 2i;
        for k = 1 : step_num_1
            A_prebor(4, sum(deg_bor(1:(k-1))+1) + 1 : sum(deg_bor(1:k)+1)) = k-1;
        end
        A_prebor(3, :) = 1;
    end
    clear tmp;
    B_prebor = zeros(num_unknCoeff_prebor, 1);
    B_prebor(1, 1) = 1;
%     A_prebor(end, :) = 1; B_prebor(end) = 0;  % Constant condition
    
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