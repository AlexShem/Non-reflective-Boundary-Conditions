function U = rodEquation_CN_wideSearch(patterns)
    par = [.5, .02, 5];
%     par = [7.8, .002, 2e4]; % Steel
%     par = [.5, h*sqrt(0.04), .05*250/4*1.60];
    D = par(2)^2;   C = par(3) * par(2)^2 / par(1);

    nu = C * tau^2 / h^4; % Параметр Куранта
    mu = D / h^2;
%     nu = .32;
%     mu = .36;
    a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*mu; s = nu/2;
    
    L = 1;  %Длина струны, м.
    T = 5; %Время, сек.
    
    Nx = single(L / h + 1); % Кол-во точек по пространству
    Nt = floor(T / tau + 1);   % Кол-во точек по времени
    
    x = linspace(0, L, Nx);
    time = linspace(0, T, Nt);
    [x, time] = meshgrid(x, time);
    true_sol = rodEquation_CN_PseudoTrueSolution([h, tau], par, [nu mu]);
    true_sol = true_sol(:, ceil(size(true_sol, 2)/2) - floor(Nx/2) : ceil(size(true_sol, 2)/2) + floor(Nx/2));
    
    u_0 = linspace(0, L, Nx)';
    u_0 = cosDistribFun_at_tao(0, 1, u_0) * 0.1;
%     u_0 = cosDistribFun_at_tao(0, 1, u_0).^4 * 0.1;
%     u_0 = (u_0/L*2*pi - pi) .* sin(u_0/L*2*pi - pi).^4 * .05;   % x * sin(x)^4 from -pi to pi 
%     u_0 = cosDistribFun_at_tao(.45, .55, u_0) * 0.1;
%     u_0 = cosDistribFun_at_tao(0, 1, u_0) .* u_0 *.25;
    %u_0 = CauchyDistribFun_at_tao_new(4.9, 5.1, u_0) * 0.1;
%     u_0 = arrayfun(@(x) 0.05, u_0);
    
    u_tau = u_0;
%     u_tau = linspace(0, L, Nx)';
%     u_tau = cosDistribFun_at_tao(0, 1, u_tau) * 0.1;
%     u_tau = cosDistribFun_at_tao(0, 1, u_tau).^4 * 0.1;
%     u_tau = cosDistribFun_at_tao(0.45, .55, u_tau) * 0.1;
%     u_tau = CauchyDistribFun_at_tao_new(0, L, u_tau) * 0.1;
%     u_tau = CauchyDistribFun_at_tao_new(0, L, u_tau) .* u_tau *.25;
%     u_tau = arrayfun(@(x) 0.05, u_tau);
%     u_tau = true_sol(2, :);
    
    U = zeros(Nt, Nx);
    U(1, :) = u_0';
    U(2, :) = u_tau';
%     ---------------------------------------------------
%     polynomDEG_bor = [0, 0, 2, 2];
%     polynomDEG_prebor = [0, 0, 2, 2];
    polynomDEG_bor = [0, 0, 10, 10];
    polynomDEG_prebor = [0, 0, 10, 10];
%     
%     polynomDEG_bor = [5,5,4,4];
%     polynomDEG_prebor = [5,5,4,4];
    
%     polynomDEG_bor = [0, 0, 30, 30];
%     polynomDEG_prebor = [0, 0, 30, 30];

% %     Good 
%     polynomDEG_bor = [4, 4, 10, 10];
%     polynomDEG_prebor = [4, 4, 10, 10];

%     polynomDEG_bor = [4, 4, 8, 8];
%     polynomDEG_prebor = [4, 4, 8, 8];

%     polynomDEG_bor = [4, 4, 28, 28];
%     polynomDEG_prebor = [4, 4, 28, 28];

%     polynomDEG_bor = patterns{1}; polynomDEG_prebor = patterns{2};
%     ---------------------------------------------------
    step_space_bor = length(polynomDEG_bor);
    step_space_prebor = length(polynomDEG_prebor);
    
    [borderDEG, preborderDEG] = PadeCoef_poly(nu, mu, polynomDEG_bor, polynomDEG_prebor);
%     [borderDEG, preborderDEG] = PadeCoef_poly_newLam(nu, mu, polynomDEG_bor, polynomDEG_prebor);
    
% load nu_01_mu_05_0088.mat;

    coefsm_bor = zeros (max(polynomDEG_bor)+1,step_space_bor);
    for i = 1:step_space_bor
        coefsm_bor (1 : (polynomDEG_bor(i) + 1), i) = 1;
    end
    coefsm_bor(coefsm_bor == 1) = borderDEG;
    
    coefsm_prebor = zeros (max(polynomDEG_prebor)+1,step_space_prebor);
    for i = 1:step_space_prebor
        coefsm_prebor (1 : (polynomDEG_prebor(i) + 1), i) = 1;
    end
    coefsm_prebor(coefsm_prebor == 1) = preborderDEG;
    
    % Создаем начальные матрицы, не учитывающие граничные коэффициенты
    U_now = diag(a*ones(1, Nx)) + diag(b * ones(1, Nx-1), 1) + diag(b * ones(1, Nx-1), -1) + diag(s * ones(1, Nx-2), 2) + diag(s * ones(1, Nx-2), -2);
    U_now(1, :) = 0; U_now(2, :) = 0;
    U_now(end, :) = 0; U_now(end-1, :) = 0;
    
    U_at_pTime = diag(d * ones(1, Nx)) + diag(g * ones(1, Nx-1), 1) + diag(g * ones(1, Nx-1), -1);
    U_at_pTime(1, :) = 0; U_at_pTime(2, :) = 0;
    U_at_pTime(end, :) = 0; U_at_pTime(end-1, :) = 0;
    
    U_at_ppTime = diag(a*ones(1, Nx)) + diag(b * ones(1, Nx-1), 1) + diag(b * ones(1, Nx-1), -1) + diag(s * ones(1, Nx-2), 2) + diag(s * ones(1, Nx-2), -2);
    U_at_ppTime(1, :) = 0; U_at_ppTime(2, :) = 0;
    U_at_ppTime(end, :) = 0; U_at_ppTime(end-1, :) = 0;
    
    % Заполняем нужные элементы матриц
    for i = 1 : step_space_prebor
        cf = coefsm_prebor(1,i);
        U_now(2, i) = cf;
        U_now(end-1, end - i + 1) = cf;
        
        cf = coefsm_prebor(2,i);
        U_at_pTime(2, i) = cf;
        U_at_pTime(end-1, end - i + 1) = cf;
        
        cf = coefsm_prebor(3,i);
        U_at_ppTime(2, i) = cf;
        U_at_ppTime(end-1, end - i + 1) = cf;
    end
    for i = 1 : step_space_bor
        cf = coefsm_bor(1,i);
        U_now(1, i) = cf;
        U_now(end, end - i + 1) = cf;
        
        cf = coefsm_bor(2,i);
        U_at_pTime(1, i) = cf;
        U_at_pTime(end, end - i + 1) = cf;
        
        cf = coefsm_bor(3,i);
        U_at_ppTime(1, i) = cf;
        U_at_ppTime(end, end - i + 1) = cf;
    end

    U_at_pTime = sparse(U_at_pTime);
    U_at_ppTime = sparse(U_at_ppTime);
    U_now = sparse(U_now);
    if max(abs(condest(U_now))) > 1e9
        error('Matrix U_now is degenerate');
    end
    
    for n = 3 : Nt
        rightPart = -U_at_pTime * U(n - 1, :)' - U_at_ppTime * U(n - 2, :)';
        addition_left_bor = 0;
        addition_left_prebor = 0;
        addition_right_bor = 0;
        addition_right_prebor = 0;
% Border
        for i = 1 : step_space_bor
            for k = 4 : polynomDEG_bor(i) + 1
                if k <= polynomDEG_bor(i)+1
                    cf = coefsm_bor(k, i);
                else 
                    break
                end
                try
                    addition_left_bor = addition_left_bor - cf * U(n - k + 1, i);
                    addition_right_bor = addition_right_bor - cf * U(n - k + 1, end - i + 1);
                catch
                    break;
                end
            end
        end
        
%PreBorder
        for i = 1 : step_space_prebor
            for k = 4 : polynomDEG_prebor(i) + 1
                if k <= polynomDEG_prebor(i)+1
                    cf = coefsm_prebor(k, i);
                else 
                    break
                end
                try
                    addition_left_prebor = addition_left_prebor - cf * U(n - k + 1, i);
                    addition_right_prebor = addition_right_prebor - cf * U(n - k + 1, end - i + 1);
                catch
                    break;
                end
            end
        end
        rightPart(1) = rightPart(1) + addition_left_bor;
        rightPart(end) = rightPart(end) + addition_right_bor;
        rightPart(2) = rightPart(2) + addition_left_prebor;
        rightPart(end-1) = rightPart(end-1) + addition_right_prebor;
%         U(n, 2:end-1) = U_now \ rightPart;
        U(n, :) = U_now \ rightPart;
        
        
        
%         U_part_left = U(n :-1: n-size(coefsm_bor, 1)+1, 1:step_space_bor);
%         border_mat = coefsm_bor .* U_part_left;
%         border_val = border_mat(1,1);
%         border_val = -sum(sum(border_mat)) + border_val;
%         U(n, 1) = border_val;
%         
%         U_part_right = U(n :-1: n-size(coefsm_bor, 1)+1, end:-1:end-step_space_bor+1);
%         border_mat = coefsm_bor .* U_part_right;
%         border_val = border_mat(1,1);
%         border_val = -sum(sum(border_mat)) + border_val;
%         U(n, end) = border_val;
        
%         plot(x(1,:), U(n,:)); axis([0, 1, 0, 0.5]);
    end

%profile report

%     U_centr = U(:, round(size(U, 2)/2));
%     figure; hold on;
%     plot(U_centr);

    patern = '';
    prepatern = '';
    patern_str = '';
    prepatern_str = '';
    for i = 1:length(polynomDEG_bor)
        patern = strcat(patern, num2str(polynomDEG_bor(i)), '_');
        patern_str = strcat(patern_str, num2str(polynomDEG_bor(i)), ', ');
    end
    patern_str = patern_str(1 : end - 1);
    for i = 1:length(polynomDEG_prebor)
        prepatern = strcat(prepatern, num2str(polynomDEG_prebor(i)), '_');
        prepatern_str = strcat(prepatern_str, num2str(polynomDEG_prebor(i)), ', ');
    end
    prepatern_str = prepatern_str(1 : end - 1);
    
%     x = linspace(0, L, Nx);
%     time = linspace(0, T, Nt);
%     [x, time] = meshgrid(x, time);
%     true_sol = wave_true_sol(x, time, c) * .1;
    
%     I = int_sum_new(U, h, tau, c);

% % Plot log_10 Wave energy
%     figure(1); hold on;
%     pp = plot(linspace(0, T, length(I)), log10(I(:,1))); xl = xlabel('t, сек.'); yl = ylabel('log_{10}(I)');
%     ttl = title({'log_{10}(Wave Energy)', ['\nu = ', num2str(nu), ', h = ', num2str(h), ', c = ', num2str(c)], ['Шаблон: ', patern_str]});
%     xl.FontSize = 18; yl.FontSize = 18;ttl.FontSize = 16; pp.LineWidth = 1.4;

% % Plot log_10 |WaveEnergy - WaveEnergy_true|
    figure(2); hold on;
    I_minus_true = int_sum_abs2(U - true_sol, h);
    pp = plot(linspace(0, T, length(I_minus_true)-1), log10(abs(I_minus_true(1:end-1)))); xl = xlabel('t, сек.'); yl = ylabel('log_{10}I_d(t)');
%     ttl = title({'log_{10}I_d(t)', ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)],...
%         ['Bor: ', patern_str, ', PreBor: ', prepatern_str]});
    ttl = title({'log_{10} I(t)', ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)], patern_str});
    axis([0 T -10 -.5]); xl.FontSize = 18; yl.FontSize = 18; ttl.FontSize = 16; pp.LineWidth = 1.4;
% %     ln = line([0.25 0.25], [-9.5 -2]); ln.LineStyle = '--'; ln.Color = 'k'; ln.LineWidth = 1.2;

% % Plot max log_10 |U - U_true|
    figure(3); hold on;
    Max_diff_U_true = max(abs(U - true_sol), [], 2);
    pp = plot(linspace(0, T, length(Max_diff_U_true)), log10(Max_diff_U_true)); xl = xlabel('t, сек.'); yl = ylabel('log_{10}M_d(t)');
%     ttl = title({'log_{10}M_d(t)', ['\nu = ', num2str(nu) '; h = ', num2str(h), '; c = ', num2str(c)], ['Шаблон: ', patern_str]});
%     ttl = title({'log_{10} max_x |\hat{u} - u| (t)', ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)],...
%         ['Bor: ', patern_str, ', PreBor: ', prepatern_str]});
    ttl = title({'log_{10} M(t)', ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)], patern_str});
    axis([0 T -6.5 -.5]); xl.FontSize = 18; yl.FontSize = 18; ttl.FontSize = 16; pp.LineWidth = 1.4;
% %     ln = line([0.25 0.25], [-6 -2.5]); ln.LineStyle = '--'; ln.Color = 'k'; ln.LineWidth = 1.2;

    patern_left = '';
    patern_left_str = '';
    for i = 1:length(polynomDEG_bor)
        patern_left = strcat(patern_left, num2str(polynomDEG_bor(i)), '_');
        patern_left_str = strcat(patern_left_str, num2str(polynomDEG_bor(i)), ',');
    end

% Lazy animation
figure(1)
for i = 1 : 5 : Nt
    subplot(1,2,1);
    plot(0:h:L, U(i,:));
    axis([0 1 -.1 .1]);
    title(num2str((i-1)*tau));
    subplot(1,2,2);
    plot(0:h:L, log10(abs(U(i,:) - true_sol(i,:))));
    title(num2str((i-1)*tau));
    axis([0 1 -6 -1])
    pause(.1);
end
    
% Animate in programm
%     figure(4);
%     for i = 1 : ceil(.005 / tau) : Nt
%         if i == 1
%             fg = plot(0:h:L, U(i, :));
%             axis([0 L -.12 .12]);
%             title({['Волна в t = ', num2str((i - 1) * tau), 'с.'],...
%                 ['Шаблон: ', patern_left_str(1 : end - 1)]});
%             xlabel('x');
%         else
%             set(fg,'YData', U(i,:));
%             set(gca, 'YLim',[-0.12 0.12]);
%             title({['Волна в t = ', num2str((i - 1) * tau), 'с.'],...
%                 ['Шаблон: ', patern_left_str(1 : end - 1)]});
%             drawnow;
%             pause(.05);
%         end
%     end
    
% %    Making animation
%     filename = strcat('abs_error_anim_', patern, 'nu_', num2str(nu), 'mu_', num2str(mu), '_t_', num2str(T), '_(no w = 1).gif');
%     figure(5);
%     for i = 1 : ceil(.0005 / tau) : Nt
% %         plot(0:h:L, U(i, :));
%         plot(0:h:L, log10(abs(U(i, :) - true_sol(i,:)))); title({'log_{10} |U_{true} - U|', ['t = ' num2str(tau * (i-1))]}); axis([0 L -16 -3]);
% %         axis([0 L -0.12 0.12]);
% %         ttl=title({['Волна в t = ', num2str((i - 1) * tau), 'с.'],...
% %             ['Шаблон: ', patern_str]}); ttl.FontSize = 16;
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