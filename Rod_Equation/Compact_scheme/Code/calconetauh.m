function [bc_exist,gd_max_H,gd_max_L2,gd_max_Ch] = calconetauh(tau, h, polynomDEG_bor, rho, R, E)
% res = '';
warning off
gd_max_H = 0;
gd_max_L2 = 0;
gd_max_Ch = 0;
try
    Nt = 5e4; % Number of time steps to be performed
    T = tau * Nt; % Integration time
    U = rodEquation_CN_htau_search(h, tau, polynomDEG_bor, T);
%     U_0 = U(1, :);
    
    bc_exist = true;
    
    %                 Energy test
    H = int_sum_derivative(U, h, tau, rho, R, E);
    H_crit = H(1);
    if all(H(ceil(.1*length(H)) : end) < H_crit)
        gd_max_H = 1;
    end
    
    % L_2 norm
    L2 = trapz(-.5 : h : .5, U.^2, 2);
    L2 = sqrt(L2);
    L2_crit = L2(1);
    if max(L2) <= L2_crit
        gd_max_L2 = 1;
    end
    
    % Chebyshev norm
    Ch = max(abs(U), [], 2);
    Ch_crit = Ch(1);
    if max(Ch) <= Ch_crit
        gd_max_Ch = 1;
    end
catch ex
    bc_exist = false;
end
%             disp(['h = ' num2str(h) ', tau = ' num2str(tau) ' done', res]);