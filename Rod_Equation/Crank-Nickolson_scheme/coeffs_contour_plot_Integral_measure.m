function [good_pattern, nu, mu] = coeffs_contour_plot_Integral_measure(deg_bor, deg_prebor)
    expr_bor = '';
    expr_prebor = '';
    for i = 1 : length(deg_bor); expr_bor = [expr_bor, num2str(deg_bor(i)), '_']; end
    for i = 1 : length(deg_prebor); expr_prebor = [expr_prebor, num2str(deg_prebor(i)), '_']; end
    
    expr = ['bor_', expr_bor(1:end-1), '_prebor_', expr_prebor(1:end-1), '.mat'];
    
    try
        load(expr);
    catch
        warning('off', 'all');
        [nu, mu] = meshgrid(linspace(1e-4, 1.5, 101), linspace(1e-4, 1.5, 101));

    end
    
    for i = 1 : size(nu, 1)
        for j = 1 : size(mu, 2)
            try
                [polyCoeff_bor, polyCoeff_prebor] = PadeCoef_poly(nu(i, j), mu(i, j), deg_bor, deg_prebor);
                coef_bor_set{i, j} = polyCoeff_bor;
                coef_prebor_set{i, j} = polyCoeff_prebor;
            catch
                warning(['No pattern for nu = ' num2str(nu(i, j)) ', mu = ' num2str(mu(i, j))]);
                coef_bor_set{i, j} = NaN;
                coef_prebor_set{i, j} = NaN;
            end
        end
    end
        
end