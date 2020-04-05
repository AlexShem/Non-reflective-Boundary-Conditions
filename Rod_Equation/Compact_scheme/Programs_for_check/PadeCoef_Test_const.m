function [polyCoeff_bor, polyCoeff_prebor] = PadeCoef_Test_const(nu, mu, deg_bor, deg_prebor)
    step_num_1 = length(deg_bor);
    step_num_2 = length(deg_prebor);
    totalDeg_1 = sum(deg_bor + 1);
    totalDeg_2 = sum(deg_prebor + 1);
    syms w;
    
    [~, ~, lam_ser_3, lam_ser_4] = lambda_series_sym(w, nu, mu,...
        max([totalDeg_1, totalDeg_2]) + 3, max([totalDeg_1, totalDeg_2]) + 3);
    
    unknCoeff_bor = {};
    for ind = 1:totalDeg_1
        unknCoeff_bor{ind} = sym(strcat('a', num2str(ind)));
    end
    unknCoeff_prebor = {};
    for ind = 1:totalDeg_2
        unknCoeff_prebor{ind} = sym(strcat('b', num2str(ind)));
    end
    
    Y_bor = struct();
    j = 1;
    eq_bor_3 = 0;
    eq_bor_4 = 0;
    for ind = 1:step_num_1
        Y_bor.(strcat('st', num2str(ind))) = '0';
        for iind = 0:deg_bor(ind)
            Y_bor.(strcat('st', num2str(ind))) = Y_bor.(strcat('st', num2str(ind))) + unknCoeff_bor{j} * w^iind;
            j = j + 1;
        end
        eq_bor_3 = eq_bor_3 + Y_bor.(strcat('st', '0' + ind)) * lam_ser_3^(ind - 1);
        eq_bor_4 = eq_bor_4 + Y_bor.(strcat('st', '0' + ind)) * lam_ser_4^(ind - 1);
    end
    
    Y_prebor = struct();
    j = 1;
    eq_prebor_3 = 0;
    eq_prebor_4 = 0;
    for ind = 1:step_num_2
        Y_prebor.(strcat('st', num2str(ind))) = '0';
        for iind = 0:deg_prebor(ind)
            Y_prebor.(strcat('st', num2str(ind))) = Y_prebor.(strcat('st', num2str(ind))) + unknCoeff_prebor{j} * w^iind;
            j = j + 1;
        end
        eq_prebor_3 = eq_prebor_3 + Y_prebor.(strcat('st', '0' + ind)) * lam_ser_3^(ind - 1);
        eq_prebor_4 = eq_prebor_4 + Y_prebor.(strcat('st', '0' + ind)) * lam_ser_4^(ind - 1);
    end
    
    
    if nu == 20100 && mu == 10000
        if isequal(deg_bor, [0,0,2,2])
            load('border_cfs.mat');
        elseif isequal(deg_bor, [1,1,2,2])
            load('border_cfs_1122.mat');
        else
            cf_bor_3 = coeffs(expand(eq_bor_3), w);
            cf_bor_4 = coeffs(expand(eq_bor_4 * w^(step_num_1 - 1)), w);
        end
    else
        cf_bor_3 = coeffs(expand(eq_bor_3), w);
        cf_bor_4 = coeffs(expand(eq_bor_4 * w^(step_num_1 - 1)), w);
    end
    if length(cf_bor_3) >= totalDeg_1
        cf_bor_3(totalDeg_1 - 1 : length(cf_bor_3)) = [];
    end
    if length(cf_bor_4) >= totalDeg_1
        cf_bor_4(totalDeg_1 - 1 : length(cf_bor_4)) = [];
    end
    
    if nu == 20100 && mu == 10000
        if isequal(deg_prebor, [0,0,2,2])
            load('preborder_cfs.mat');
        elseif isequal(deg_bor, [1,1,2,2])
            load('preborder_cfs_1122.mat');
        else
            cf_prebor_3 = coeffs(expand(eq_prebor_3), w);
            cf_prebor_4 = coeffs(expand(eq_prebor_4 * w^(step_num_2 - 1)), w); 
        end
    else
        cf_prebor_3 = coeffs(expand(eq_prebor_3), w);
        cf_prebor_4 = coeffs(expand(eq_prebor_4 * w^(step_num_2 - 1)), w); 
    end
    if length(cf_prebor_3) >= totalDeg_2
        cf_prebor_3(totalDeg_2 - 1: length(cf_prebor_3)) = [];
    end
    if length(cf_prebor_4) >= totalDeg_2
        cf_prebor_4(totalDeg_2 - 1 : length(cf_prebor_4)) = [];
    end
    
    A_bor = zeros(length(unknCoeff_bor));
    B_bor = zeros(length(unknCoeff_bor), 1);
    B_bor(1, 1) = 1;
%     B_bor(2, 1) = 1;
    for i = 1 : length(unknCoeff_bor)
        if i <= totalDeg_1 / 2
            cf_copy = cf_bor_3;
        else
            cf_copy = cf_bor_4;
        end
        subVar = zeros(1, length(unknCoeff_bor));
        subVar(i) = 1;
        %cf_copy = double(subs(cf_copy, unknCoeff, subVar));
%         A_bor(2:end, i) = double(subs(cf_copy, unknCoeff_bor, subVar))';
        A_bor(3:end, i) = double(subs(cf_copy, unknCoeff_bor, subVar))';
    end
    A_bor(1, 1) = 1;
    A_bor(2, deg_bor(1) + 2) = 1;
    
    A_prebor = zeros(length(unknCoeff_prebor));
    B_prebor = zeros(length(unknCoeff_prebor), 1);
    B_prebor(1, 1) = 0;
    B_prebor(2, 1) = 1;
    for i = 1:length(unknCoeff_prebor)
        if i <= totalDeg_2 / 2
            cf_copy = cf_prebor_3;
        else
            cf_copy = cf_prebor_4;
        end
        subVar = zeros(1, length(unknCoeff_prebor));
        subVar(i) = 1;
        %cf_copy = double(subs(cf_copy, unknCoeff, subVar));
%         A_prebor(2:end, i) = double(subs(cf_copy, unknCoeff_prebor, subVar))';
        A_prebor(3:end, i) = double(subs(cf_copy, unknCoeff_prebor, subVar))';
    end
    A_prebor(1, 1) = 1;
    A_prebor(2, deg_prebor(1) + 2) = 1;
    
%     deg_count = 1;
%     for i = 1 : step_num
% % Use lambda_series(w, nu, mu) intead
%         A(end, deg_count : deg_count + degrees(i)) = lam_1_true_fun(nu, 1)^(step_num - i);
%         deg_count = deg_count + degrees(i) + 1;
%     end
    
    eig_val_bor = eig(A_bor); eig_bor_min = min(abs(eig_val_bor));
    eig_val_prebor = eig(A_prebor); eig_prebor_min = min(abs(eig_val_prebor));

    if eig_bor_min <= 1e-5 && eig_prebor_min <= 1e-5
        error('Both sets of border points have no Pade approximation');
    elseif eig_bor_min <= 1e-5
        error('Border point set has no Pade aproximation');
    elseif eig_prebor_min <= 1e-5
        error('Preborder point set has no Pade aproximation');
    end
    
    polyCoeff_bor = linsolve(A_bor, B_bor);
    polyCoeff_prebor = linsolve(A_prebor, B_prebor);
end