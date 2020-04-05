function [polyCoeff_bor, polyCoeff_prebor] = PadeCoef_Test_const(nu, mu, deg_bor, deg_prebor)
    step_num_1 = length(deg_bor);
    step_num_2 = length(deg_prebor);
    totalDeg_1 = sum(deg_bor + 1);
    totalDeg_2 = sum(deg_prebor + 1);
    maxTotalDeg = max([totalDeg_1, totalDeg_2]);
    syms w;
%     assume(w, 'real');
    
    [~, ~, lam_ser_3, lam_ser_4] = lambda_series_sym(w, nu, mu, maxTotalDeg, maxTotalDeg);
    
    lam_3 = sym2poly(lam_ser_3);    display('lam_3 parsed to poly');
    lam_4 = sym2poly(lam_ser_4);    display('lam_4 parsed to poly');
    lam_3_real = poly2sym(lam_3(end-maxTotalDeg : end), w);
    lam_4_real = poly2sym(lam_4(end-maxTotalDeg : end), w);

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
        eq_bor_3 = eq_bor_3 + Y_bor.(strcat('st', '0' + ind)) * lam_3_real^(ind - 1);
        eq_bor_4 = eq_bor_4 + Y_bor.(strcat('st', '0' + ind)) * lam_4_real^(ind - 1);
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
        eq_prebor_3 = eq_prebor_3 + Y_prebor.(strcat('st', '0' + ind)) * lam_3_real^(ind - 1);
        eq_prebor_4 = eq_prebor_4 + Y_prebor.(strcat('st', '0' + ind)) * lam_4_real^(ind - 1);
    end
    
    bor_1_coef_3 = coeffs(collect(eq_bor_3, w), w); cf_bor_1 = bor_1_coef_3; display('1. border lam 1 real ready');
    bor_1_coef_4 = coeffs(collect(eq_bor_4, w), w); cf_bor_2 = bor_1_coef_4; display('2. border lam 1 imag ready');
    
    if length(cf_bor_1) >= totalDeg_1
        cf_bor_1(totalDeg_1 - 1 : length(cf_bor_1)) = [];
    end
    if length(cf_bor_2) >= totalDeg_1
        cf_bor_2(totalDeg_1 - 1 : length(cf_bor_2)) = [];
    end
    
    prebor_1_coef_3 = coeffs(collect(eq_prebor_3, w), w); cf_prebor_1 = prebor_1_coef_3; display('3. preborder lam 1 real ready');
    prebor_1_coef_4 = coeffs(collect(eq_prebor_4, w), w); cf_prebor_2 = prebor_1_coef_4; display('4. preborder lam 1 imag ready');
    
    if length(cf_prebor_1) >= totalDeg_2
        cf_prebor_1(totalDeg_2: length(cf_prebor_1)) = [];
    end
    if length(cf_prebor_2) >= totalDeg_2
        cf_prebor_2(totalDeg_2 : length(cf_prebor_2)) = [];
    end
    
    A_bor = zeros(length(unknCoeff_bor));
    B_bor = zeros(length(unknCoeff_bor), 1);
    a_bor_size = size(A_bor, 1);
    B_bor(1, 1) = 1;
%     B_bor(2, 1) = 1;
    for i = 1 : length(unknCoeff_bor)
        cf_copy_1 = cf_bor_1;
        cf_copy_2 = cf_bor_2;
        subVar = zeros(1, length(unknCoeff_bor));
        subVar(i) = 1;
        cfs_1 = double(subs(cf_copy_1, unknCoeff_bor, subVar))';
        cfs_2 = double(subs(cf_copy_2, unknCoeff_bor, subVar))';
        
        A_bor(3 : 3-1+ceil((a_bor_size-2)/2), i) = cfs_1(1 : ceil((a_bor_size-2)/2));
        A_bor(3+ceil((a_bor_size-2)/2) : end, i) = cfs_2(1 : floor((a_bor_size-2)/2));
    end
    A_bor(1, 1) = 1;
    A_bor(2, deg_bor(1) + 2) = 1;
    
    A_prebor = zeros(length(unknCoeff_prebor));
    B_prebor = zeros(length(unknCoeff_prebor), 1);
    a_prebor_size = size(A_prebor, 1);
    B_prebor(1, 1) = 0;
    B_prebor(2, 1) = 1;
    for i = 1:length(unknCoeff_prebor)
        cf_copy_1 = cf_prebor_1;
        cf_copy_2 = cf_prebor_2;
        subVar = zeros(1, length(unknCoeff_prebor));
        subVar(i) = 1;
        
        cfs_1 = double(subs(cf_copy_1, unknCoeff_prebor, subVar))';
        cfs_2 = double(subs(cf_copy_2, unknCoeff_prebor, subVar))';

        A_prebor(3 : 3-1+ceil((a_prebor_size-2)/2), i) = cfs_1(1 : ceil((a_prebor_size-2)/2));
        A_prebor(3+ceil((a_prebor_size-2)/2) : end, i) = cfs_2(1 : floor((a_prebor_size-2)/2));
    end
    A_prebor(1, 1) = 1;
    A_prebor(2,:) = 0; A_prebor(2, deg_prebor(1) + 2) = 1;
    
% % Условие на константе
%     A_bor(end, :) = 1;
%     A_prebor(end, :) = 1;
    
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