function U = rodEquation_CN(patterns)
%profile on;
    h = 0.005; tau = 0.00001;

%     Set the parameters
%   par = [rho, R, E]
    par = [7860, 1e-3, 210e9];  % %For Steel 
    D = par(2)^2;   C = par(3) * par(2)^2 / par(1);

    nu = C * tau^2 / h^4; % Параметр Куранта
    mu = D / h^2;
    a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*mu; s = nu/2;
    
    L = 1;  %Длина струны, м.
    T = .015; %Время, сек.
    
    Nx = single(L / h + 1); % Кол-во точек по пространству
    Nt = floor(T / tau + 1);   % Кол-во точек по времени
    
    x = linspace(0, L, Nx);
    time = linspace(0, T, Nt);
    [x, time] = meshgrid(x, time);
    true_sol = rodEquation_CN_PseudoTrueSolution([h, tau], par, [nu mu], T);
    true_sol = true_sol(:, ceil(size(true_sol, 2)/2) - floor(Nx/2) : ceil(size(true_sol, 2)/2) + floor(Nx/2));
    
    u_0 = linspace(0, L, Nx)';
    u_0 = cosDistribFun_at_tao(0, 1, u_0);
    
    u_tau = u_0;
    
    U = zeros(Nt, Nx);
    U(1, :) = u_0';
    U(2, :) = u_tau';
%     ---------------------------------------------------
%     polynomDEG_bor = [5, 5, 18, 18];
%     polynomDEG_prebor = [5, 5, 18, 18];

%     polynomDEG_bor = [5, 5, 14, 14];
%     polynomDEG_prebor = [5, 5, 14, 14];
    
%     polynomDEG_bor = [5, 5, 10, 10];
%     polynomDEG_prebor = [5, 5, 10, 10];
    
%     polynomDEG_bor = [4, 4, 8, 8];
%     polynomDEG_prebor = [4, 4, 8, 8];
    
%     polynomDEG_bor = [3, 3, 9, 9];
%     polynomDEG_prebor = [3, 3, 9, 9];
    
%     polynomDEG_bor = [6, 4, 12, 8];
%     polynomDEG_prebor = [6, 4, 12, 8];
    
%    %Philipp work
%     polynomDEG_bor = [4, 4, 8, 8];
%     polynomDEG_bor_step = [5, 5, 10, 10];
%     polynomDEG_bor = [3, 5, 8, 8];    
%     polynomDEG_bor_step = [5, 3, 8, 8];
%     polynomDEG_bor_step = [4, 4, 8, 8];
%     polynomDEG_bor_step = [5, 5, 11, 13];
%     polynomDEG_bor = [5, 5, 11, 13];
%     polynomDEG_bor = [5, 7, 13, 17]; %%%%%%%%%%%
%     polynomDEG_bor = [5, 7, 13, 21];
%     polynomDEG_bor = [4, 4, 12, 12];
    polynomDEG_bor = [3, 5, 9, 9];

%   % Test
%     polynomDEG_bor = [4, 4, 4, 4];

    polynomDEG_prebor = polynomDEG_bor;
%     ---------------------------------------------------
    step_space_bor = length(polynomDEG_bor);
    step_space_prebor = length(polynomDEG_prebor);
    
    [borderDEG, preborderDEG] = PadeCoef_poly(nu, mu, polynomDEG_bor, polynomDEG_prebor);
    
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
%         error('Matrix U_now is degenerate');
%     end

    nStepT_bor = max(polynomDEG_bor) + 1;
    nStepT_prebor = max(polynomDEG_prebor) + 1;
    opersRight = cell(max(nStepT_bor, nStepT_prebor) - 1, 1);
    opersRight{1} = U_at_pTime; opersRight{2} = U_at_ppTime;
    tstep = max(nStepT_bor, nStepT_prebor);
    
    
%     % Main Cycle
    for n = 3 : Nt
        if n <= tstep && n > 3
            operT = zeros(Nx);
            operT([1 2], 1:step_space_bor) = [coefsm_bor(n, :); coefsm_prebor(n, :)];
            operT([end end-1], end:-1:end-step_space_bor+1) = [coefsm_bor(n, :); coefsm_prebor(n, :)];
            opersRight{n-1} = sparse(operT);
        end
        rightPart = 0;
        for k = 1 : min(n-1, tstep-1)
            rightPart = rightPart + opersRight{k} * U(n-k, :).';
        end
        U(n, :) = -U_now \ rightPart;
    end

%profile report

    patern_left = '';
    patern_left_str = '';
    for i = 1:length(polynomDEG_bor)
        patern_left = strcat(patern_left, num2str(polynomDEG_bor(i)), '_');
        patern_left_str = strcat(patern_left_str, num2str(polynomDEG_bor(i)), ',');
    end
    patern_left_str = patern_left_str(1:end-1);
    patern_right = '';
    patern_right_str = '';
    for i = 1:length(polynomDEG_prebor)
        patern_right = strcat(patern_right, num2str(polynomDEG_prebor(i)), '_');
        patern_right_str = strcat(patern_right_str, num2str(polynomDEG_prebor(i)), ',');
    end
    patern_right_str = patern_right_str(1:end-1);

% %     % Analyse the operator ----------------------------------------------
%     operSize = Nx * length(opersRight);
%     E = eye(Nx);
%     lower = eye(length(opersRight)-1);
%     lower(:, end+1) = 0;
%     lower = kron(lower, E);
%     uper = zeros(Nx, operSize);
%     for k = 1 : length(opersRight)
%         uper(:, Nx*(k-1)+1 : Nx*k) = opersRight{k};
%     end
%     oper = [uper; lower];
% %     lower = eye(length(opersRight) - 1);
% %     operLeft = kron(lower, eye(Nx));
% %     operLeft = [full(U_now), zeros(Nx, operSize-Nx); zeros(operSize - Nx, Nx), operLeft];
%     oper = [full(U_now), zeros(Nx, operSize-Nx); zeros(operSize-Nx, Nx), eye(operSize - Nx)] \ oper;
%     [V, lam] = eig(oper); lam = diag(lam);
%     figure(5)
%     plot(sort(abs(lam)), '.b');
%     xlabel('Номер собственного числа');
%     title({'Собственные значения оператора перехода',...
%         ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h)],...
%         [patern_left_str, '; ', patern_right_str]})
%     
%     % Plotting eigenvectors
%     m = max(abs(lam));
%     ind = abs(lam) == m;
%     vComp = V(:, ind);
%     if size(vComp, 2) == 2
%         v = [sum(vComp)/2, diff(vComp)/2i];
%     else
%         v = vComp;
%     end
%     figure(7);
%     for k = 1 : size(v, 2)
%         plot(v(:, k)); hold on;
%     end
%     hold off;
%     xlabel('Номер координаты');
%     title(['Координаты собственного вектора для |\lambda| = ', num2str(abs(lam(ind)))]);
%     
%     
%     [~, sind] = sort(abs(lam));
%     slam = lam(sind);
%     vComp = V(:, sind);
%     vComp = vComp(:, end:-1:end-1);
% %     if size(vComp, 2) == 2
% %         v = [sum(vComp)/2, diff(vComp)/2i];
% %     else
% %         v = vComp;
% %     end
%     v = vComp;
%     figure(8);
%     for k = 1 : size(v, 2)
%         plot(v(:, k)); hold on;
%     end
%     hold off;
%     xlabel('Номер координаты');
%     title(['Координаты собственного вектора для |\lambda| = ', num2str(abs(lam(ind)))]);
% %     %--------------------------------------------------------------------

% % Plot log_10 |WaveEnergy - WaveEnergy_true|
    figure(1); hold on;
    I_minus_true = int_sum_abs2(U - true_sol, h);
    pp = plot(linspace(0, T, length(I_minus_true)-1), log10(abs(I_minus_true(1:end-1)))); xl = xlabel('t, сек.'); yl = ylabel('log_{10}I_d(t)');
%     ttl = title({'log_{10}I_d(t)', ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)],...
%         ['Bor: ', patern_str, ', PreBor: ', prepatern_str]});
    ttl = title({'log_{10} I(t)', ...
        ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)], ...
        [patern_left_str, ' ; ', patern_right_str]});
    axis([0 T -4 0]); xl.FontSize = 18; yl.FontSize = 18; ttl.FontSize = 16; pp.LineWidth = 1.4;
% %     ln = line([0.25 0.25], [-9.5 -2]); ln.LineStyle = '--'; ln.Color = 'k'; ln.LineWidth = 1.2;
    grid on;

% % Plot max log_10 |U - U_true|
    figure(2); hold on;
    Max_diff_U_true = max(abs(U - true_sol), [], 2);
    pp = plot(linspace(0, T, length(Max_diff_U_true)), log10(Max_diff_U_true)); xl = xlabel('t, сек.'); yl = ylabel('log_{10}M_d(t)');
%     ttl = title({'log_{10}M_d(t)', ['\nu = ', num2str(nu) '; h = ', num2str(h), '; c = ', num2str(c)], ['Шаблон: ', patern_str]});
%     ttl = title({'log_{10} max_x |\hat{u} - u| (t)', ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)],...
%         ['Bor: ', patern_str, ', PreBor: ', prepatern_str]});
    ttl = title({'log_{10} M(t)', ...
        ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)], ...
        [patern_left_str, ' ; ', patern_right_str]});
    axis([0 T -4 0]); xl.FontSize = 18; yl.FontSize = 18; ttl.FontSize = 16; pp.LineWidth = 1.4;
% %     ln = line([0.25 0.25], [-6 -2.5]); ln.LineStyle = '--'; ln.Color = 'k'; ln.LineWidth = 1.2;
    grid on;

% Lazy animation
figure(6)
for i = 1 : 1 : Nt
    subplot(1,2,1);
    plot(0:h:L, U(i,:));
    axis([0 1 -.15 1]);
    title(num2str((i-1)*tau));
    subplot(1,2,2);
    plot(0:h:L, log10(abs(U(i,:) - true_sol(i,:))));
    title(num2str((i-1)*tau));
    axis([0 1 -6 0])
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
%     filename = strcat('abs_error_anim_', patern_left_str, 'nu_', num2str(nu), 'mu_', num2str(mu), '_t_', num2str(T), '.gif');
%     figure(5);
%     for i = 1 : ceil(.0005 / tau) : Nt
% %         plot(0:h:L, U(i, :));
%         plot(0:h:L, log10(abs(U(i, :) - true_sol(i,:)))); title({'log_{10} |U_{true} - U|', ['t = ' num2str(tau * (i-1))]}); axis([0 L -6 0]);
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