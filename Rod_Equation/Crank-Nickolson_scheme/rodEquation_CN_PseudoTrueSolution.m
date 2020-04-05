function U = rodEquation_CN_PseudoTrueSolution(steps, par, kurant, T)
%     load Compact_Pade_Coeffs;
%profile on
% INITIAL h = 0.005; tau = 0.001;

%     h = 0.005; tau = 0.001;   %Precise
%     h = 0.00025; tau = 0.00025;
%     h = 0.001; tau = 0.0002;

%     Set the parameters
%   par = [rho, R, E]
%     par = [1, sqrt(1.5)*h/10, 1];  % nu = 0.015; mu = 0.015;
%     par = [1, sqrt(1.5/3)*h, 100/6/(10/3)];  % nu = 0.1; mu = 0.5;
%     par = [1, sqrt(1.5/3)*h, 100/6/(1/3)];  % nu = 1; mu = 0.5;
%     par = [.5, .001, 5];
%     D = par(2)^2;   C = par(3) * par(2)^2 / par(1);
%     D=0;
%     C=0;

    h = steps(1); tau = steps(2);
    D = par(2)^2;   C = par(3) * par(2)^2 / par(1);

%     nu = C * tau^2 / h^4; % Параметр Куранта
%     mu = D / h^2;

    nu = kurant(1); mu = kurant(2);
    a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*mu; s = nu/2;
    
    L = 10;  %Длина струны, м.
%     L = 20;
%     T = .05; %Время, сек.
%     T = 20;
%     T = 10;
    
    Nx = single(L / h + 1); % Кол-во точек по пространству
    Nt = floor(T / tau + 1);   % Кол-во точек по времени
%     Nt = floor(T / tau);
    
%     x = linspace(0, L, Nx);
%     time = linspace(0, T, Nt);
%     [x, time] = meshgrid(x, time);
%     true_sol = wave_true_sol(x, time, D) * .1;
    
    u_0 = linspace(0, L, Nx)';
    u_0 = linspace(-L/2, L/2, Nx)'; u_0 = cosDistribFun_at_tao(-.5, .5, u_0) .* u_0;
%     u_0 = cosDistribFun_at_tao(9.5, 10.5, u_0);
%     u_0 = cosDistribFun_at_tao(10 - pi/2, 10 + pi/2, u_0).^2;
%     u_0 = cosDistribFun_at_tao(4.5, 5.5, u_0).^4 * 0.1;
%     u_0 = linspace(-L/2, L/2, Nx);
%     dom = abs(u_0) <= pi/2;
%     u_0(~dom) = 0;
%     u_0(dom) = cos(u_0(dom)).^2;
%     for i = 1 : length(u_0) % x * sin(x)^4 *0.05 from -pi to pi 
%         if u_0(i) >= 4.5 && u_0(i) <= 5.5
%             u_0(i) = (u_0(i) - 5)*2*pi * sin((u_0(i) - 5)*2*pi)^4 * .05;
%         else
%             u_0(i) = 0;
%         end
%     end
%     u_0 = cosDistribFun_at_tao(4.95, 5.05, u_0) * 0.1;
%     u_0 = cosDistribFun_at_tao(0, 1, u_0) .* u_0 *.25;
%     u_0 = CauchyDistribFun_at_tao_new(4.9, 5.1, u_0) * 0.1;
%     u_0 = arrayfun(@(x) 0.05, u_0);
    
    u_tau = u_0;
%     u_tau = linspace(0, L, Nx)';
%     u_tau = cosDistribFun_at_tao(4.5, 5.5, u_tau) * 0.1;
%     u_tau = cosDistribFun_at_tao(4.5, 5.5, u_tau).^4 * 0.1;
%     u_tau = CauchyDistribFun_at_tao_new(4.95, 5.05, u_tau) * 0.1;
%     u_tau = CauchyDistribFun_at_tao_new(0, L, u_tau) .* u_tau *.25;
%     u_tau = arrayfun(@(x) 0.05, u_tau);
%     u_tau = true_sol(2, :);
    
    U = zeros(Nt, Nx);
    U(1, :) = u_0';
    U(2, :) = u_tau';
%     ---------------------------------------------------
%     polynomDEG_bor = [0, 0, 2, 2];
%     polynomDEG_prebor = [0, 0, 2, 2];
%     polynomDEG_bor = [0, 0, 3, 3];
%     polynomDEG_prebor = [0, 0, 3, 3];

    polynomDEG_bor = [0, 0, 2, 2];
    polynomDEG_prebor = [0, 0, 2, 2];

%     polynomDEG_bor = patterns{1}; polynomDEG_prebor = patterns{2};
%     ---------------------------------------------------
    step_space_bor = length(polynomDEG_bor);
    step_space_prebor = length(polynomDEG_prebor);
    
%     [borderDEG, preborderDEG] = PadeCoef_Test_const(nu, mu, polynomDEG_bor, polynomDEG_prebor);
%     [borderDEG, preborderDEG] = PadeCoef_Test_complex(nu, mu, polynomDEG_bor, polynomDEG_prebor);
    
% load nu_01_mu_05_0088.mat;

%     coefsm_bor = zeros (max(polynomDEG_bor)+1,step_space_bor);
%     for i = 1:step_space_bor
%         coefsm_bor (1 : (polynomDEG_bor(i) + 1), i) = 1;
%     end
%     coefsm_bor(coefsm_bor == 1) = borderDEG;
% % TMP-------------------
% % coefsm_bor(1,2) = 0;
% % coefsm_bor(:,3:4)=-1;
%     
%     coefsm_prebor = zeros (max(polynomDEG_prebor)+1,step_space_prebor);
%     for i = 1:step_space_prebor
%         coefsm_prebor (1 : (polynomDEG_prebor(i) + 1), i) = 1;
%     end
%     coefsm_prebor(coefsm_prebor == 1) = preborderDEG;
% % TMP-------------------
% % coefsm_prebor(1,1) = 0;
% % coefsm_prebor(:,3:4)=-1;
    
    % Создаем начальные матрицы, не учитывающие граничные коэффициенты
%     U_now = diag(ones(1, Nx-2)) + diag(a * ones(1, Nx-3), 1) + diag(a * ones(1, Nx-3), -1);
%     U_now(1, :) = 0;
%     U_now(end, :) = 0;
    U_now = diag(a*ones(1, Nx)) + diag(b * ones(1, Nx-1), 1) + diag(b * ones(1, Nx-1), -1) + diag(s * ones(1, Nx-2), 2) + diag(s * ones(1, Nx-2), -2);
%     U_now(1, :) = 0; U_now(2, :) = 0;
%     U_now(end, :) = 0; U_now(end-1, :) = 0;
%  Diriclet Boundary
    U_now([1,2,end-1, end],:) = []; U_now(:,[1,2,end-1,end]) = [];

    U_at_pTime = diag(d * ones(1, Nx)) + diag(g * ones(1, Nx-1), 1) + diag(g * ones(1, Nx-1), -1);
%     U_at_pTime(1, :) = 0; U_at_pTime(2, :) = 0;
%     U_at_pTime(end, :) = 0; U_at_pTime(end-1, :) = 0;
    U_at_pTime([1,2,end-1, end],:) = []; U_at_pTime(:,[1,2,end-1,end]) = [];
    
    U_at_ppTime = diag(a*ones(1, Nx)) + diag(b * ones(1, Nx-1), 1) + diag(b * ones(1, Nx-1), -1) + diag(s * ones(1, Nx-2), 2) + diag(s * ones(1, Nx-2), -2);
%     U_at_ppTime(1, :) = 0; U_at_ppTime(2, :) = 0;
%     U_at_ppTime(end, :) = 0; U_at_ppTime(end-1, :) = 0;
    U_at_ppTime([1,2,end-1, end],:) = []; U_at_ppTime(:,[1,2,end-1,end]) = [];

    % Заполняем нужные элементы матриц
%     for i = 1 : step_space_prebor
%         cf = coefsm_prebor(1,i);
%         if i > 1
%             U_now(1, i-1) = cf;
%             U_now(end, end - i + 2) = cf;
%         end
%         
%         cf = coefsm_prebor(2,i);
%         U_at_pTime(1, i) = cf;
%         U_at_pTime(end, end - i + 1) = cf;
%         
%         cf = coefsm_prebor(3,i);
%         U_at_ppTime(1, i) = cf;
%         U_at_ppTime(end, end - i + 1) = cf;
%     end

%     for i = 1 : step_space_prebor
%         cf = coefsm_prebor(1,i);
%         U_now(2, i) = cf;
%         U_now(end-1, end - i + 1) = cf;
%         
%         cf = coefsm_prebor(2,i);
%         U_at_pTime(2, i) = cf;
%         U_at_pTime(end-1, end - i + 1) = cf;
%         
%         cf = coefsm_prebor(3,i);
%         U_at_ppTime(2, i) = cf;
%         U_at_ppTime(end-1, end - i + 1) = cf;
%     end
%     for i = 1 : step_space_bor
%         cf = coefsm_bor(1,i);
%         U_now(1, i) = cf;
%         U_now(end, end - i + 1) = cf;
%         
%         cf = coefsm_bor(2,i);
%         U_at_pTime(1, i) = cf;
%         U_at_pTime(end, end - i + 1) = cf;
%         
%         cf = coefsm_bor(3,i);
%         U_at_ppTime(1, i) = cf;
%         U_at_ppTime(end, end - i + 1) = cf;
%     end

%     U_at_pTime = sparse(U_at_pTime);
%     U_at_ppTime = sparse(U_at_ppTime);
%     U_now = sparse(U_now);
%     if max(abs(condest(U_now))) > 1e9
%         error('Matrix U_now is degenerate');
%     end
    
    U_now = sparse(U_now);
    U_at_pTime = sparse(U_at_pTime);
    U_at_ppTime = sparse(U_at_ppTime);
    for n = 3 : Nt
        rightPart = -U_at_pTime * U(n - 1, 3:end-2)' - U_at_ppTime * U(n - 2, 3:end-2)';
        U(n, 3:end-2) = U_now \ rightPart;
    end

%profile report

%     U_centr = U(:, round(size(U, 2)/2));
%     figure; hold on;
%     plot(U_centr);

    patern = '';
    patern_str = '';
    for i = 1:length(polynomDEG_bor)
        patern = strcat(patern, num2str(polynomDEG_bor(i)), '_');
        patern_str = strcat(patern_str, num2str(polynomDEG_bor(i)), ', ');
    end
    patern_str = patern_str(1 : end - 1);
    
%     x = linspace(0, L, Nx);
%     time = linspace(0, T, Nt);
%     [x, time] = meshgrid(x, time);
%     true_sol = wave_true_sol(x, time, c) * .1;
    
%     I = int_sum_new(U, h, tau, c);

% Plot log_10 Wave energy
%     figure(1); hold on;
%     pp = plot(linspace(0, T, length(I)), log10(I(:,1))); xl = xlabel('t, сек.'); yl = ylabel('log_{10}(I)');
%     ttl = title({'log_{10}(Wave Energy)', ['\nu = ', num2str(nu), ', h = ', num2str(h), ', c = ', num2str(c)], ['Шаблон: ', patern_str]});
%     xl.FontSize = 18; yl.FontSize = 18;ttl.FontSize = 16; pp.LineWidth = 1.4;

% Plot log_10 |WaveEnergy - WaveEnergy_true|
%     figure(2); hold on;
%     I_minus_true = int_sum_new(U - true_sol, h, tau, c);
%     pp = plot(linspace(0, T, length(I)-1), log10(abs(I_minus_true(1:end-1)))); xl = xlabel('t, сек.'); yl = ylabel('log_{10}I_d(t)');
%     ttl = title({'log_{10}I_d(t)', ['\nu = ', num2str(nu) '; h = ', num2str(h), '; c = ', num2str(c)], ['Шаблон: ', patern_str]});
%     axis([0 T -10 -4]); xl.FontSize = 18; yl.FontSize = 18; ttl.FontSize = 16; pp.LineWidth = 1.4;
% %     ln = line([0.25 0.25], [-9.5 -2]); ln.LineStyle = '--'; ln.Color = 'k'; ln.LineWidth = 1.2;

% Plot max log_10 |U - U_true|
%     figure(3); hold on;
%     Max_diff_U_true = max(abs(U - true_sol), [], 2);
%     pp = plot(linspace(0, T, length(I)), log10(Max_diff_U_true)); xl = xlabel('t, сек.'); yl = ylabel('log_{10}M_d(t)');
%     ttl = title({'log_{10}M_d(t)', ['\nu = ', num2str(nu) '; h = ', num2str(h), '; c = ', num2str(c)], ['Шаблон: ', patern_str]});
%     axis([0 T -6.5 -3]); xl.FontSize = 18; yl.FontSize = 18; ttl.FontSize = 16; pp.LineWidth = 1.4;
% %     ln = line([0.25 0.25], [-6 -2.5]); ln.LineStyle = '--'; ln.Color = 'k'; ln.LineWidth = 1.2;

    patern_left = '';
    patern_left_str = '';
    for i = 1:length(polynomDEG_bor)
        patern_left = strcat(patern_left, num2str(polynomDEG_bor(i)), '_');
        patern_left_str = strcat(patern_left_str, num2str(polynomDEG_bor(i)), ',');
    end

% Animate in programm
%     figure(4);
%     for i = 1 : ceil(.05 / tau) : Nt
%         if i == 1
%             fg = plot(0:h:L, U(i, :));
% %             axis([4.5 5.5 -.12 .12]);
%             axis([0 10 -.12 1]);
%             title({['Волна в t = ', num2str((i - 1) * tau), 'с.'],...
%                 ['Шаблон: ', patern_left_str(1 : end - 1)]});
%             xlabel('x');
%         else
%             set(fg,'YData', U(i,:));
%             set(gca, 'YLim',[-0.12 1]);
%             title({['Волна в t = ', num2str((i - 1) * tau), 'с.'],...
%                 ['Шаблон: ', patern_left_str(1 : end - 1)]});
%             drawnow;
%             pause(.05);
%         end
%     end

% %    Making animation
%     filename = strcat('anim_', patern, 'nu_', num2str(nu), 'mu_', num2str(mu), '_t_', num2str(T), '_(no w = 1).gif');
%     figure(5);
%     for i = 1 : ceil(.0005 / tau) : Nt
%         plot(0:h:L, U(i, :));
%         axis([0 L -0.12 0.12]);
%         ttl=title({['Волна в t = ', num2str((i - 1) * tau), 'с.'],...
%             ['Шаблон: ', patern_str]}); ttl.FontSize = 16;
% %         ttl = title({'Волна в ', ['t = ', num2str((i - 1) * tau), 'с.']}); ttl.FontSize = 16;
%         xlabel('x', 'FontSize', 16);
%         drawnow;
%         frame = getframe(5);
%         im = frame2im(frame);
%         [A, map] = rgb2ind(im, 256);
%         if i == 1
%             imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.001);
%         else
%             imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0);
%         end
%     end
end