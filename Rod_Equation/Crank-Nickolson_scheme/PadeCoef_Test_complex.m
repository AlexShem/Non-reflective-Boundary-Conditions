function [polyCoeff_bor, polyCoeff_prebor] = PadeCoef_Test_complex(nu, mu, deg_bor, deg_prebor)
    step_num_1 = length(deg_bor);
    step_num_2 = length(deg_prebor);
    totalDeg_1 = sum(deg_bor + 1);
    totalDeg_2 = sum(deg_prebor + 1);
    maxTotalDeg = max([totalDeg_1, totalDeg_2]);
    syms w;
    assume(w, 'real');
    
    syms a b;
    
    [~, ~, lam_ser_3, ~] = lambda_series_sym(w, nu, mu, maxTotalDeg, maxTotalDeg);
    
%     lam_1_real = real(lam_ser_1);
%     lam_1_imag = imag(lam_ser_1);
    
    lam_3 = sym2poly(lam_ser_3);
    lam_ser_3_real = poly2sym(real(lam_3(end-maxTotalDeg : end)), w);
    lam_ser_3_imag = poly2sym(imag(lam_3(end-maxTotalDeg : end)), w);

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
        
        eq_for_lam_plus = simplify(((a+1i*b)^(ind-1) + (a-1i*b)^(ind-1)) / 2);
        if ind == 1
            eq_for_lam_minus = sym(1);
        else
            eq_for_lam_minus = simplify(((a+1i*b)^(ind-1) - (a-1i*b)^(ind-1)) / (2i));
        end
        eq_bor_3 = eq_bor_3 + Y_bor.(strcat('st', '0' + ind)) * subs(eq_for_lam_plus, [a b], [lam_ser_3_real, lam_ser_3_imag]);
        eq_bor_4 = eq_bor_4 + Y_bor.(strcat('st', '0' + ind)) * subs(eq_for_lam_minus, [a b], [lam_ser_3_real, lam_ser_3_imag]);
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
        
        eq_for_lam_plus = simplify(((a+1i*b)^(ind-1) + (a-1i*b)^(ind-1)) / 2);
        if ind == 1
            eq_for_lam_minus = sym(1);
        else
            eq_for_lam_minus = simplify(((a+1i*b)^(ind-1) - (a-1i*b)^(ind-1)) / (2i));
        end
        eq_prebor_3 = eq_prebor_3 + Y_prebor.(strcat('st', '0' + ind)) * subs(eq_for_lam_plus, [a,b], [lam_ser_3_real, lam_ser_3_imag]);
        eq_prebor_4 = eq_prebor_4 + Y_prebor.(strcat('st', '0' + ind)) * subs(eq_for_lam_minus, [a,b], [lam_ser_3_real, lam_ser_3_imag]);
    end
    
%     eq_bor = eq_bor_3 + eq_bor_4;
%     bor_coef = coeffs(collect(eq_bor, w), w); display('Border Ready.');
%     if length(bor_coef) >= totalDeg_1
%         bor_coef(totalDeg_1 - 1 : end) = [];
%     end
    
    bor_3_coef = coeffs(collect(eq_bor_3, w), w);
    bor_4_coef = coeffs(collect(eq_bor_4, w), w);
    cf_bor_3 = bor_3_coef; display('1. border lam 3 ready');
    cf_bor_4 = bor_4_coef; display('2. border lam 4 ready');
    if length(cf_bor_3) >= totalDeg_1
        cf_bor_3(totalDeg_1 - 1 : length(cf_bor_3)) = [];
    end
    if length(cf_bor_4) >= totalDeg_1
        cf_bor_4(totalDeg_1 - 1 : length(cf_bor_4)) = [];
    end
    
%     eq_prebor = eq_prebor_3 + eq_prebor_4;
%     prebor_coef = coeffs(collect(eq_prebor, w), w); display('Preborder Ready.');
%     if length(prebor_coef) >= totalDeg_1
%         prebor_coef(totalDeg_2 : end) = [];
%     end

    prebor_3_coef = coeffs(collect(eq_prebor_3, w), w);
    prebor_4_coef = coeffs(collect(eq_prebor_4, w), w);
    cf_prebor_3 = prebor_3_coef; display('3. preborder lam 3 ready');
    cf_prebor_4 = prebor_4_coef; display('4. preborder lam 4 ready');
    if length(cf_prebor_3) >= totalDeg_2
        cf_prebor_3(totalDeg_2: length(cf_prebor_3)) = [];
    end
    if length(cf_prebor_4) >= totalDeg_2
        cf_prebor_4(totalDeg_2 : length(cf_prebor_4)) = [];
    end
    
    A_bor = zeros(length(unknCoeff_bor));
    a_bor_size = size(A_bor, 1);
    B_bor = zeros(length(unknCoeff_bor), 1);
    B_bor(1, 1) = 1;
%     B_bor(2, 1) = 1;
    for i = 1 : length(unknCoeff_bor)
%         if i <= totalDeg_1 / 2
%             cf_copy = cf_bor_3;
%             subVar = zeros(1, length(unknCoeff_bor));
%             subVar(i) = 1;
%         else
%             cf_copy = cf_bor_4;
%             subVar = zeros(1, length(unknCoeff_bor));
%             subVar(i) = 1;
%         end
        cf_copy_1 = cf_bor_3;
        cf_copy_2 = cf_bor_4;
        subVar = zeros(1, length(unknCoeff_bor));
        subVar(i) = 1;
        %cf_copy = double(subs(cf_copy, unknCoeff, subVar));
%         A_bor(2:end, i) = double(subs(cf_copy, unknCoeff_bor, subVar))';
%         A_bor(3:end, i) = double(subs(cf_copy, unknCoeff_bor, subVar))';
        
        cfs_1 = double(subs(cf_copy_1, unknCoeff_bor, subVar))';
        cfs_2 = double(subs(cf_copy_2, unknCoeff_bor, subVar))';

        A_bor(3 : 3-1+ceil((a_bor_size-2)/2), i) = cfs_1(1 : ceil((a_bor_size-2)/2));
        A_bor(3+ceil((a_bor_size-2)/2) : end, i) = cfs_2(1 : floor((a_bor_size-2)/2));
    end
    A_bor(1, 1) = 1;
    A_bor(2, deg_bor(1) + 2) = 1;
    
    A_prebor = zeros(length(unknCoeff_prebor));
    a_prebor_size = size(A_prebor, 1);
    B_prebor = zeros(length(unknCoeff_prebor), 1);
    B_prebor(1, 1) = 0;
    B_prebor(2, 1) = 1;
    for i = 1:length(unknCoeff_prebor)
%         if i <= totalDeg_2 / 2
%             cf_copy = cf_prebor_3;
%             subVar = zeros(1, length(unknCoeff_prebor));
%             subVar(i) = 1;
%         else
%             cf_copy = cf_prebor_4;
%             subVar = zeros(1, length(unknCoeff_prebor));
%             subVar(i - ceil(totalDeg_2/2)) = 1;
%         end
%         cf_copy = prebor_coef;
        cf_copy_1 = cf_prebor_3;
        cf_copy_2 = cf_prebor_4;
        subVar = zeros(1, length(unknCoeff_prebor));
        subVar(i) = 1;
        %cf_copy = double(subs(cf_copy, unknCoeff, subVar));
%         A_prebor(2:end, i) = double(subs(cf_copy, unknCoeff_prebor, subVar)).';
%         A_prebor(3:end, i) = double(subs(cf_copy, unknCoeff_prebor, subVar)).';
%         A_prebor(2:end, i) = double(subs(cf_copy, unknCoeff_prebor, subVar)).';

        cfs_1 = double(subs(cf_copy_1, unknCoeff_prebor, subVar))';
        cfs_2 = double(subs(cf_copy_2, unknCoeff_prebor, subVar))';

        A_prebor(3 : 3-1+ceil((a_prebor_size-2)/2), i) = cfs_1(1 : ceil((a_prebor_size-2)/2));
        A_prebor(3+ceil((a_prebor_size-2)/2) : end, i) = cfs_2(1 : floor((a_prebor_size-2)/2));
    end
    A_prebor(1, 1) = 1;
    A_prebor(2,:) = 0; A_prebor(2, deg_prebor(1) + 2) = 1;
    
%     deg_count = 1;
%     for i = 1 : step_num
% Use lambda_series(w, nu, mu) intead
%         A(end, deg_count : deg_count + degrees(i)) = lam_1_true_fun(nu, 1)^(step_num - i);
%         deg_count = deg_count + degrees(i) + 1;
%     end

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