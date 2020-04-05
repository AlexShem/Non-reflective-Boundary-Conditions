function U = rodEquation_Compact_Dirichlet()
%     load Compact_Pade_Coeffs;
%profile on
    
    h = 0.005; tau = 0.001;   %Precise
%     h = 0.00025; tau = 0.00025;
%     h = 0.001; tau = 0.0002;

%     Set the parameters
%   par = [rho, R, E]
%     par = [1, 1, 100];
%     par = [1, 100*h, 201*h^2/(100*tau^2)];
    par = [1, 100*h, 9*h^2/(10*tau^2)];
    D = par(2)^2;   C = par(3) * par(2)^2 / par(1);

    nu = C * tau^2 / h^4; % Параметр Куранта
    mu = D / h^2;
    a = 3/(12*mu + 4) - .5; b = 9*nu/(3*mu + 1) - 2; c = 1 - (12*nu + 3) / (6*mu + 2); d = 3*nu / (6*mu + 2); 
    
    L = 1;  %Длина струны, м.
    T = 1; %Время, сек.
    
    Nx = single(L / h + 1); % Кол-во точек по пространству
    Nt = single(T / tau + 1);   % Кол-во точек по времени
    
    x = linspace(0, L, Nx);
    time = linspace(0, T, Nt);
    [x, time] = meshgrid(x, time);
    true_sol = wave_true_sol(x, time, D) * .1;
    
    u_0 = linspace(0, L, Nx)';
    u_0 = true_sol(1,:);
%     u_0 = cosDistribFun_at_tao(0, 1, u_0) * 0.1;
    %u_0 = CauchyDistribFun_at_tao_new(0, L, u_0) * 0.1;
%     u_0 = CauchyDistribFun_at_tao_new(0.2, 0.8, u_0) * 0.1;
    %u_0 = CauchyDistribFun_at_tao_new(4.9, 5.1, u_0) * 0.1;
    %u_0 = arrayfun(@(x) 0.05, u_0);
    
%     u_tau = u_0;
%     u_tau = linspace(0, L, Nx)';
%     u_tau = cosDistribFun_at_tao(0, 1, u_tau) * 0.1;
    %u_tau = CauchyDistribFun_at_tao_new(0, L, u_tau) * 0.1;
%     u_tau = CauchyDistribFun_at_tao_new(0.2, 0.8, u_tau) * 0.1;
    %u_tau = CauchyDistribFun_at_tao_new(4.9, 5.1, u_tau) * 0.1;
    %u_tau = arrayfun(@(x) 0.05, u_tao);
    u_tau = true_sol(2, :);
    
    U = zeros(Nt, Nx);
    U(1, :) = u_0';
    U(2, :) = u_tau';
    
    % Fro Dirichlet
%     U_now = diag(ones(1, Nx-4)) + diag(a * ones(1, Nx-5), 1) + diag(a * ones(1, Nx-5), -1);
%     
%     U_at_pTime = diag(b * ones(1, Nx - 4)) + diag(c * ones(1, Nx-5), 1) + ...
%         diag(c * ones(1, Nx-5), -1) + diag(d * ones(1, Nx-6), +2) + diag(d * ones(1, Nx-6), -2);
%     
%     U_at_ppTime = diag(ones(1, Nx-4)) + diag(a * ones(1, Nx-5), 1) + diag(a * ones(1, Nx-5), -1);
    
    % For Neuman
    U_now = diag(ones(1, Nx-4)) + diag(a * ones(1, Nx-5), 1) + diag(a * ones(1, Nx-5), -1);
    
    U_at_pTime = diag(b * ones(1, Nx - 4)) + diag(c * ones(1, Nx-5), 1) + ...
        diag(c * ones(1, Nx-5), -1) + diag(d * ones(1, Nx-6), +2) + diag(d * ones(1, Nx-6), -2);
    
    U_at_ppTime = diag(ones(1, Nx-4)) + diag(a * ones(1, Nx-5), 1) + diag(a * ones(1, Nx-5), -1);
        
%     U_at_pTime = sparse(U_at_pTime);
%     U_at_ppTime = sparse(U_at_ppTime);
%     U_now = sparse(U_now);
%     if max(abs(condest(U_now))) > 1e9
%         error('Matrix U_now is degenerate');
%     end
    
    for n = 3 : Nt
        % For Diriclet
        rightPart = -U_at_pTime * U(n - 1, 3:end-2)' - U_at_ppTime * U(n - 2, 3:end-2)';
        U(n, 3:end-2) = U_now \ rightPart;
        
        % For Neuman
%         rightPart = -U_at_pTime * U(n - 1, 3:end-2)' - U_at_ppTime * U(n - 2, 3:end-2)';
%         U(n, 3:end-2) = U_now \ rightPart;
%         U(n, 2) = U(n, 3); U(n, 1) = U(n, 3);
%         U(n, end-1) = U(n, end-2); U(n, end) = U(n, end-2);
    end

%profile report

%     U_centr = U(:, round(size(U, 2)/2));
%     figure; hold on;
%     plot(U_centr);
        
    I = int_sum_new(U, h, tau, c);

% Plot log_10 Wave energy
%     figure(1); hold on;
%     pp = plot(linspace(0, T, length(I)), log10(I(:,1))); xl = xlabel('t, сек.'); yl = ylabel('log_{10}(I)');
%     ttl = title({'log_{10}(Wave Energy)', ['\nu = ', num2str(nu), ', h = ', num2str(h), ', c = ', num2str(c)], ['Шаблон: ', patern_str]});
%     xl.FontSize = 18; yl.FontSize = 18;ttl.FontSize = 16; pp.LineWidth = 1.4;

% Plot log_10 |WaveEnergy - WaveEnergy_true|
    figure(2); hold on;
    I_minus_true = int_sum_new(U - true_sol, h, tau, c);
    pp = plot(linspace(0, T, length(I)-1), log10(abs(I_minus_true(1:end-1)))); xl = xlabel('t, сек.'); yl = ylabel('log_{10}I_d(t)');
    ttl = title({'log_{10}I_d(t)', ['\nu = ', num2str(nu) '; h = ', num2str(h), '; c = ', num2str(c)], ['Шаблон: ', patern_str]});
    axis([0 T -10 -4]); xl.FontSize = 18; yl.FontSize = 18; ttl.FontSize = 16; pp.LineWidth = 1.4;
% %     ln = line([0.25 0.25], [-9.5 -2]); ln.LineStyle = '--'; ln.Color = 'k'; ln.LineWidth = 1.2;

% Plot max log_10 |U - U_true|
    figure(3); hold on;
    Max_diff_U_true = max(abs(U - true_sol), [], 2);
    pp = plot(linspace(0, T, length(I)), log10(Max_diff_U_true)); xl = xlabel('t, сек.'); yl = ylabel('log_{10}M_d(t)');
    ttl = title({'log_{10}M_d(t)', ['\nu = ', num2str(nu) '; h = ', num2str(h), '; c = ', num2str(c)], ['Шаблон: ', patern_str]});
    axis([0 T -6.5 -3]); xl.FontSize = 18; yl.FontSize = 18; ttl.FontSize = 16; pp.LineWidth = 1.4;
% %     ln = line([0.25 0.25], [-6 -2.5]); ln.LineStyle = '--'; ln.Color = 'k'; ln.LineWidth = 1.2;

% Animate in programm
    figure(4);
    for i = 1 : ceil(.005 / tau) : Nt
        if i == 1
            fg = plot(0:h:L, U(i, :));
            axis([0 L -0.12 0.12]);
            title({['Волна в t = ', num2str((i - 1) * tau), 'с.'],...
                ['Шаблон: ', patern_left_str(1 : end - 1)]});
            xlabel('x');
        else
            set(fg,'YData', U(i,:));
            set(gca, 'YLim',[-0.12 0.12]);
            title({['Волна в t = ', num2str((i - 1) * tau), 'с.'],...
                ['Шаблон: ', patern_left_str(1 : end - 1)]});
            drawnow;
            pause(0.01);
        end
    end
end