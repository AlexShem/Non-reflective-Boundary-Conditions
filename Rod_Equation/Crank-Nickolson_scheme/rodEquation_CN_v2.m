function U = rodEquation_CN_v2(patterns)
%profile on;
    % % Best [4 5 8 8], const condition; h = 0.02; tau = 1.6e-4; Steel
    % parameters
    
%     h = 0.05469; tau = 5.017e-05;
    h = .02; tau = 1.6e-4;
%     h = .02/2; tau = 1.6e-4/4;

    h = .01; tau = 1.6e-4;

%     Set the parameters
%   par = [rho, R, E]
    par = [7860, 1e-3, 210e9];  % %For Steel 
%     par = [7860, 1.1e-3, 210e9];  % %Slight change
    D = par(2)^2;   C = par(3) * par(2)^2 / par(1);

    nu = C * tau^2 / h^4; % Параметр Куранта
    mu = D / h^2;
    a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*mu; s = nu/2;
    
    L = 1;  %Длина струны, м.
%     T = .001; %Время, сек.
    T = .3; %Время, сек.
%     T = .01*10; %Время, сек.
    
%     Nx = floor(L / h + 1); % Кол-во точек по пространству
%     Nt = floor(T / tau + 1);   % Кол-во точек по времени
    
    Nx = length(-L/2 : h : L/2); % Кол-во точек по пространству
    Nt = length(0 : tau : T);   % Кол-во точек по времени
            
    u_0 = linspace(0, L, Nx)';
    x = -L/2 : h : L/2;
%     Symmetric initial condition
    [u_0, sigma] = GaussianDistrib(-.3, .3, x); 
    u_0 = u_0 .* x; u_0 = u_0(:);

%     Shifted initial condition
%     [u_0, sigma] = GaussianDistrib(-.3, .3, x-.1); 
%     u_0 = u_0 .* (x - .1); u_0 = u_0(:);
    
%     u_0 = .125*(cos(pi*(x-.1)/.1)+1) .* (x >= 0 & x <= 0.2) - ...
%         .125*(cos(pi*(x-.1)/.1)+1) .* (x < 0 & x >= -0.2); u_0 = u_0(:);
%     u_0 = .125*(cos(pi*(x)/.1)+1) .* (x >= 0.1 & x <= 0.3) - ...
%         .125*(cos(pi*(x)/.1)+1) .* (x < 0.1 & x >= -0.1); u_0 = u_0(:);

    h_min = min(u_0);
    h_max = max(u_0);
    
    u_tau = u_0;
%     f = GaussianDistrib(-.3, .3, (-20*L) : h : (20*L)); f = f .* ((-20*L) : h : (20*L)); f = f(:);
%     u_tau = u_tau_function(h, D, C, f); len_tau = length(u_tau);
%     u_tau = u_tau(ceil((len_tau+1)/2 - length(u_0)/2) : floor((len_tau+1)/2 + length(u_0)/2));
%     u_tau = u_0 + u_tau * tau^2 / 2; u_tau = u_tau(:);
    
    U = zeros(Nt, Nx);
    U(1, :) = u_0';
    U(2, :) = u_tau';
%     ---------------------------------------------------
%   % NO CONST
%     polynomDEG_bor = [4, 4, 8, 8]; constant_condition = 0;
%     polynomDEG_bor = [8, 8, 6, 6]; constant_condition = 0;
    polynomDEG_bor = [7, 6, 6, 9]; constant_condition = 0;
% %     polynomDEG_bor = [8, 8, 8, 8]; constant_condition = 0;
%     polynomDEG_bor = [5, 3, 8, 8]; constant_condition = 0;
%     polynomDEG_bor = [5, 3, 9, 7]; constant_condition = 0;  
    
%     polynomDEG_bor = [4, 5, 8, 8]; constant_condition = 1;
%     polynomDEG_bor = [5, 3, 10, 7]; constant_condition = 1;   
    
%   % CONST
%     polynomDEG_bor = [4, 5, 8, 8]; constant_condition = 1;
%     polynomDEG_bor = [6, 8, 9, 10]; constant_condition = 1;
%     polynomDEG_bor = [6, 9, 7, 7]; constant_condition = 1;
%     polynomDEG_bor = [7, 6, 6, 8]; constant_condition = 1;

%     polynomDEG_bor = patterns;
%     constant_condition = 0;
    polynomDEG_prebor = polynomDEG_bor;
%     ---------------------------------------------------
    step_space_bor = length(polynomDEG_bor);
    step_space_prebor = length(polynomDEG_prebor);
    
    x_condition = 0;
    
%     [borderDEG, preborderDEG] = ...
%         PadeCoef_poly(nu, mu, polynomDEG_bor, polynomDEG_prebor, true);
    [borderDEG, preborderDEG] = ...
        PadeCoef_poly_const_x(nu, mu, polynomDEG_bor, polynomDEG_prebor, ...
        constant_condition, x_condition);
%     [borderDEG, preborderDEG] = ...
%         PadeCoef_poly_const_x_1w2(nu, mu, polynomDEG_bor, polynomDEG_prebor, ...
%         constant_condition, x_condition);
    
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
    
%     % Заполняем нужные элементы матриц
    U_now(1:2, 1:step_space_bor) = [coefsm_bor(1, :); coefsm_prebor(1, :)];
    U_now(end:-1:end-1, end:-1:end-step_space_bor+1) = [coefsm_bor(1, :); coefsm_prebor(1, :)];
    U_at_pTime(1:2, 1:step_space_bor) = [coefsm_bor(2, :); coefsm_prebor(2, :)];
    U_at_pTime(end:-1:end-1, end:-1:end-step_space_bor+1) = [coefsm_bor(2, :); coefsm_prebor(2, :)];
    U_at_ppTime(1:2, 1:step_space_bor) = [coefsm_bor(3, :); coefsm_prebor(3, :)];
    U_at_ppTime(end:-1:end-1, end:-1:end-step_space_bor+1) = [coefsm_bor(3, :); coefsm_prebor(3, :)];
    
    U_at_pTime = sparse(U_at_pTime);
    U_at_ppTime = sparse(U_at_ppTime);
    U_now = sparse(U_now);
%     if max(abs(condest(U_now))) > 1e9
%         warning('Matrix U_now is degenerate');
%     end
%     U_now_inv = inv(U_now);

    nStepT_bor = max(polynomDEG_bor) + 1;
    nStepT_prebor = max(polynomDEG_prebor) + 1;
    
%     % Main Cycle
    for n = 3 : Nt
%         if n <= tstep && n > 3
%             operT = zeros(Nx);
%             operT([1 2], 1:step_space_bor) = [coefsm_bor(n, :); coefsm_prebor(n, :)];
%             operT([end end-1], end:-1:end-step_space_bor+1) = [coefsm_bor(n, :); coefsm_prebor(n, :)];
%             opersRight{n-1} = sparse(operT);
%         end
%         rightPart = 0;
%         for k = 1 : min(n-1, tstep-1)
%             rightPart = rightPart + opersRight{k} * U(n-k, :).';
%         end
        
        % MAKER A BORDER HERE
        cf_rem_bor = coefsm_bor(4 : min(nStepT_bor, n), :);
        cf_rem_prebor = coefsm_prebor(4 : min(nStepT_prebor, n), :);
        U_prev_left_bor = U(n-3 : -1 : max(n-nStepT_bor+1, 1), 1 : step_space_bor);
        U_prev_left_prebor = U(n-3 : -1 : max(n-nStepT_prebor+1, 1), 1 : step_space_prebor);
        U_prev_right_bor = U(n-3 : -1 : max(n-nStepT_bor+1, 1), end : -1 : end-step_space_bor+1);
        U_prev_right_prebor = U(n-3 : -1 : max(n-nStepT_prebor+1, 1), end : -1 : end-step_space_prebor+1);
        add_left_bor = cf_rem_bor .* U_prev_left_bor;
        add_left_bor = sum(add_left_bor(:));
        add_left_prebor = cf_rem_prebor .* U_prev_left_prebor; add_left_prebor = sum(add_left_prebor(:));
        add_right_bor = cf_rem_bor .* U_prev_right_bor; add_right_bor = sum(add_right_bor(:));
        add_right_prebor = cf_rem_prebor .* U_prev_right_prebor; add_right_prebor = sum(add_right_prebor(:));
        
        rightPart = U_at_pTime * U(n-1, :).' + U_at_ppTime * U(n-2, :).';
        rightPart(1) = rightPart(1) + add_left_bor;
        rightPart(2) = rightPart(2) + add_left_prebor;
        rightPart(end-1) = rightPart(end-1) + add_right_prebor;
        rightPart(end) = rightPart(end) + add_right_bor;

        U(n, :) = - U_now \ rightPart;
    end

%profile report

    patern_left = '';
    patern_left_str = '';
    for i = 1:length(polynomDEG_bor)
        patern_left = strcat(patern_left, num2str(polynomDEG_bor(i)), '_');
        patern_left_str = strcat(patern_left_str, num2str(polynomDEG_bor(i)), '-');
    end
    patern_left_str = patern_left_str(1:end-1);
    patern_right = '';
    patern_right_str = '';
    for i = 1:length(polynomDEG_prebor)
        patern_right = strcat(patern_right, num2str(polynomDEG_prebor(i)), '_');
        patern_right_str = strcat(patern_right_str, num2str(polynomDEG_prebor(i)), '-');
    end
    patern_right_str = patern_right_str(1:end-1);

%     % True Solution
    x_domain = linspace(0, L, Nx);
    time = linspace(0, T, Nt);
   
    tic
    true_sol = rodEquation_CN_PseudoTrueSolution_v2([h, tau], par, [nu mu], T, u_0, u_tau, x_domain);
    toc
    
    center = size(true_sol, 2) / 2 + .5;
    LB = ceil(center - Nx/2);
    RB = floor(center + Nx/2);
    true_sol = true_sol(:, LB:RB);

% % Plot Energy norm
    f1 = figure(11);
    hold on;
    [H_int, K_int, P_int] = int_sum_derivative(U - true_sol, h, tau, par(1), par(2), par(3));
    H = sqrt(H_int);
    pp = plot(linspace(0, T, length(H)), log10(H));
    xl = xlabel('$t$ (s)', 'Interpreter', 'latex', 'FontSize', 18);
    yl = ylabel('$\log_{10} \mathcal{H} [u_j(t) - u_j^*(t)]$', 'Interpreter', 'latex', 'FontSize', 18);
    axis([0 T -1 4.5]); pp.LineWidth = 1.4;
    ax1 = f1.CurrentAxes; ax1.FontSize = 14;
    
% % % Plot C norm
    f2 = figure(12);
    hold on;
    Max_diff_U_true = max(abs(U - true_sol), [], 2);
    pp = plot(linspace(0, T, length(Max_diff_U_true)), log10(Max_diff_U_true));
    xl = xlabel('$t$ (s)', 'Interpreter', 'latex', 'FontSize', 18);
    yl = ylabel('$\log_{10} \max_j |u_j(t) - u_j^*(t)|$', 'Interpreter', 'latex', 'FontSize', 18);
    axis([0 T -5.5 -1]); pp.LineWidth = 1.4;
    ax2 = f2.CurrentAxes; ax2.FontSize = 14;
    
% % Plor L2 norm
    f3 = figure(13);
    hold on;
    L2_diff_U_true = sqrt(trapz(-L/2 : h : L/2, (U - true_sol).^2, 2));
    pp = plot(linspace(0, T, length(L2_diff_U_true)), log10(L2_diff_U_true));
    xl = xlabel('$t$ (s)', 'Interpreter', 'latex', 'FontSize', 18);
    yl = ylabel('$\log_{10} \|u_j(t) - u_j^*(t)\|_2$', 'Interpreter', 'latex', 'FontSize', 18);
    axis([0 T -5.5 -1]); pp.LineWidth = 1.4;
    ax3 = f3.CurrentAxes; ax3.FontSize = 14;

% % Lazy animation
% figure(14)
% pause(1);
% for i = 1 : 10 : Nt
%     subplot(1,2,1);
%     plot(0:h:L, U(i,:)); hold on;
%     plot(0:h:L, true_sol(i, :), '--k'); hold off;
% %     axis([0 L .01*h_min .01*h_max]);
%     axis([0 L h_min h_max]);
%     title(['t = ' num2str((i-1)*tau)]);
%     
%     subplot(1,2,2);
%     plot(0:h:L, log10(abs(U(i,:) - true_sol(i,:))));
%     title(['t = ' num2str((i-1)*tau)]);
%     axis([0 L -6 0])
%     pause(.01);
% end
    
    
% %    Making animation
%     filename = strcat('anim_', patern_left_str, '.gif');
%     figure(15);
%     pause(1);
%     for i = 1 : 10 : Nt
% 
%         subplot(1,2,1);
%         plot(-L/2 : h : L/2, U(i, :)); hold on;
%         plot(-L/2 : h : L/2, true_sol(i, :), '--k'); hold off;
%         axis([-L/2 L/2 h_min h_max]);
%         xlabel('$x\,(m)$', 'FontSize', 16, 'Interpreter', 'latex');
%         lg = legend('$u$', '$u^*$'); lg.FontSize = 18; lg.Location = 'northeast'; lg.Interpreter = 'latex';
%         title(patern_left_str);
%         
%         subplot(1,2,2);
%         plot(-L/2 : h : L/2, log10(abs(U(i, :) - true_sol(i,:))));
%         title(['t = ' num2str(tau * (i-1))]);
%         axis([-L/2 L/2 -6 0]);
%         xlabel('$x\,(m)$', 'FontSize', 16, 'Interpreter', 'latex');
% %         yl = ylabel('$log_{10} |u - u^*|$', 'Interpreter', 'latex'); yl.FontSize = 16;
%         lg = legend('$\log_{10} |u - u^*|$', 'Interpreter', 'latex'); lg.FontSize = 16; lg.Interpreter = 'latex';
% 
%         drawnow;
%         frame = getframe(15);
%         im = frame2im(frame);
%         [A, map] = rgb2ind(im, 256);
%         if i == 1
%             imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.001);
%         else
%             imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0);
%         end
%     end

end