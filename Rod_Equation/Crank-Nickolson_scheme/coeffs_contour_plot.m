function [coef_bor_set, coef_prebor_set, nu, mu] = coeffs_contour_plot(deg_bor, deg_prebor)%, coef_ind)
    
    expr_bor = '';
    expr_prebor = '';
    for i = 1 : length(deg_bor); expr_bor = [expr_bor, num2str(deg_bor(i)), '_']; end
    for i = 1 : length(deg_prebor); expr_prebor = [expr_prebor, num2str(deg_prebor(i)), '_']; end
    expr_bor = ['bor_', expr_bor(1:end-1), '.mat']; expr_prebor = ['prebor_', expr_prebor(1:end-1), '.mat'];
    
    try
        load(expr_bor);
        load(expr_prebor);
    catch    
        warning('off', 'all');
        [nu, mu] = meshgrid(linspace(1e-4, 1.5, 101), linspace(1e-4, 1.5, 101));
    %     deg_bor = [6,6,5,5];
    %     deg_prebor = [6,6,5,5];

        coef_bor_set = cell(size(nu));
        coef_prebor_set = cell(size(nu));

    %     coef_ind = 17;
    %     coef_bor = NaN(size(nu));
    %     coef_prebor = NaN(size(nu));

        for i = 1 : size(nu, 1)
            for j = 1 : size(mu, 2)
                try
                    [polyCoeff_bor, polyCoeff_prebor] = PadeCoef_poly(nu(i, j), mu(i, j), deg_bor, deg_prebor);
                    coef_bor_set{i, j} = polyCoeff_bor;
                    coef_prebor_set{i, j} = polyCoeff_prebor;
    %                 coef_bor(i, j) = polyCoeff_bor(coef_ind);
    %                 coef_prebor(i, j) = polyCoeff_prebor(coef_ind);
                catch
                    warning(['No pattern for nu = ' num2str(nu(i, j)) ', mu = ' num2str(mu(i, j))]);
                    coef_bor_set{i, j} = NaN;
                    coef_prebor_set{i, j} = NaN;
                end
            end
        end
        warning('on', 'all');
        
        save(expr_bor, 'coef_bor_set', 'nu', 'mu');
        save(expr_prebor, 'coef_prebor_set', 'nu', 'mu');
    end
    
%     figure(6)
%     contour(nu, mu, coef_bor_set, 'ShowText', 'on');
%     % contour(nu, mu, log10(abs(coef_bor)), 'ShowText', 'on');
%     xlabel('\nu'); ylabel('\mu');
%     title({['Border, coeff number ' num2str(coef_ind)], deg_bor})
end