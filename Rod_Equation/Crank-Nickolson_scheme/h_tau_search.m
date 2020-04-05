warning('off','all');

par = [7860, 1e-3, 210e9];
rho = par(1); R = par(2); E = par(3);

% % No constant condition
degrees = {[4, 4, 8, 8];
%     [8, 8, 6, 6];
%     [5, 3, 8, 8];
    [5, 3, 9, 7]};

%%---------- Precise calculations ----------
% [h, tau] = meshgrid(...
%     .0005 : .000125 : .1,...
%     .00001 : .0000004 : .00025);
[h, tau] = meshgrid(...
    linspace(.0001, .07, 750),...
    linspace(.00001, .0006, 750));

%%----------Testing the programm ----------
% [h, tau] = meshgrid(...
%     linspace(.0001, .07, 30),...
%     linspace(.00001, .0006, 30));
% [h, tau] = meshgrid(...
%     linspace(.0001, .25, 40),...
%     linspace(.00001, .0025, 40));
% [h, tau] = meshgrid(...
%     linspace(.0001, .5, 20),...
%     linspace(.00001, .005, 20));

nx = size(h, 2); nt = size(tau, 1);
gd_max_H = zeros(size(h));
gd_max_L2 = zeros(size(h));
gd_max_Ch = zeros(size(h));
bc_exist = false(size(h));

for k = 1 : length(degrees)
    tic
    polynomDEG_bor = degrees{k};
    for j = 1 : nx
        for i = 1 : nt
            res = '';
%             polynomDEG_bor = degrees{k};
            try
                Nt = 1e5; % Number of time steps to be performed
                T = tau(i, j) * Nt; % Integration time
                U = rodEquation_CN_htau_search(h(i, j), tau(i, j), polynomDEG_bor, T);
%                 U_0 = U(1, :);
                
                bc_exist(i, j) = true;
                
%                 Energy test
                H = int_sum_derivative(U, h(i, j), tau(i, j), rho, R, E);
                H_crit = H(1);
                if all(H(ceil(.1*length(H)) : end) < H_crit)
                    gd_max_H(i, j) = 1;
                end

                % L_2 norm
                L2 = trapz(-.5 : h(i, j) : .5, U.^2, 2);
                L2 = sqrt(L2);
                L2_crit = L2(1);
                if max(L2) <= L2_crit
                    gd_max_L2(i, j) = 1;
                end
                
                % Chebyshev norm
                Ch = max(abs(U), [], 2);
                Ch_crit = Ch(1);
                if max(Ch) <= Ch_crit
                    gd_max_Ch(i, j) = 1;
                end
            catch ex
                bc_exist(i, j) = false;
            end
            %             disp(['h = ' num2str(h(i,j)) ', tau = ' num2str(tau(i,j)) ' done', res]);
        end
        
        if j == floor(nx/2)
            disp(['50% for ', num2str(polynomDEG_bor)]);
        elseif j == floor(nx/4)
            disp(['25% for ', num2str(polynomDEG_bor)]);
        elseif j == floor(.75*nx)
            disp(['75% for ', num2str(polynomDEG_bor)]);
        elseif j == nx
            disp(['100% for ', num2str(polynomDEG_bor)]);
        end
    end
    toc
    
    figure(3*k - 2);
    imagesc(h(1,:), tau(:,1), gd_max_H); axis xy; colormap([0 0 0; 1 1 1]);
    xlabel('$h$ (m)', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('$\tau$ (s)', 'Interpreter', 'latex', 'FontSize', 15);
    title(['H ', num2str(polynomDEG_bor)]);
    
    figure(3*k - 1);
    imagesc(h(1,:), tau(:,1), gd_max_Ch); axis xy; colormap([0 0 0; 1 1 1]);
    xlabel('$h$ (m)', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('$\tau$ (s)', 'Interpreter', 'latex', 'FontSize', 15);
    title(['C max ', num2str(polynomDEG_bor)]);
    
    figure(3*k);
    imagesc(h(1,:), tau(:,1), gd_max_L2); axis xy; colormap([0 0 0; 1 1 1]);
    xlabel('$h$ (m)', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('$\tau$ (s)', 'Interpreter', 'latex', 'FontSize', 15);
    title(['L^2 ', num2str(polynomDEG_bor)]);
    
    filename = [num2str(polynomDEG_bor), 'no_const.mat'];
    save(filename, 'h', 'tau', 'gd_max_L2', 'gd_max_Ch', 'bc_exist');
    filename = [num2str(polynomDEG_bor), 'no_const_energy.mat'];
    save(filename, 'h', 'tau', 'gd_max_H');
end