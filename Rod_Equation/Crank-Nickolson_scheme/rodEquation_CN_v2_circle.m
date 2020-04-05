function U = rodEquation_CN_v2_circle(h, tau)
%     h = 0.005*4; tau = 0.00001*16;
    h = 0.5; tau = 0.002;

%     Set the parameters
%   par = [rho, R, E]
    par = [7860, 1e-3, 210e9];  % %For Steel 
%     par = [7860, 1.1e-3, 210e9];  % %Slight change
    D = par(2)^2;   C = par(3) * par(2)^2 / par(1);

    nu = C * tau^2 / h^4; % Параметр Куранта
    mu = D / h^2;
    a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*mu; s = nu/2;
    
    L = 1;
    T = .04*18;
    
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
    
    % Создаем начальные матрицы, не учитывающие граничные коэффициенты
%     U_now = diag(a*ones(1, Nx)) + diag(b * ones(1, Nx-1), 1) + diag(b * ones(1, Nx-1), -1) + diag(s * ones(1, Nx-2), 2) + diag(s * ones(1, Nx-2), -2);
    U_now = spdiags(ones(Nx, 1) .* [s b a b s], -2:2, Nx, Nx);
    U_now(end-1, 1) = s; U_now(2, end) = s;
    U_now(end, 1) = b; U_now(end, 2) = s;
    U_now(1, end) = b; U_now(1, end-1) = s;
    
%     U_at_pTime = diag(d * ones(1, Nx)) + diag(g * ones(1, Nx-1), 1) + diag(g * ones(1, Nx-1), -1);
    U_at_pTime = spdiags(ones(Nx, 1) .* [g d g], -1:1, Nx, Nx);
    U_at_pTime(1, end) = g;
    U_at_pTime(end, 1) = g;
    
    U_at_ppTime = U_now;
    
%     % Main Cycle
    for n = 3 : Nt
        U(n, :) = - U_now \ (U_at_pTime*U(n-1, :).' + U_at_ppTime*U(n-2, :).');
%         U(n, :) = - U_now \ U_at_pTime * U(n-1, :).'; U(n, :) = U(n, :) - U(n-2, :);
%         if n == floor(Nt / 4)
%             disp('25% is done');
%         elseif n == floor(Nt / 2)
%             disp('50% is done');
%         elseif n == floor(3 * Nt / 4)
%             disp('75% is done');
%         end
    end
    
% % Plot max log_10 |U - U_true|
%     f2 = figure(12);
%     hold on;
%     Max_diff_U_true = max(abs(U - true_sol), [], 2);
%     pp = plot(linspace(0, T, length(Max_diff_U_true)), log10(Max_diff_U_true));
%     xl = xlabel('$t$ (s)', 'Interpreter', 'latex', 'FontSize', 18);
%     yl = ylabel('$\log_{10} \max_j |u_j(t) - u_j^*(t)|$', 'Interpreter', 'latex', 'FontSize', 18);
%     axis([0 T -5.5 -1]); pp.LineWidth = 1.4;
%     ax2 = f2.CurrentAxes; ax2.FontSize = 14;
    
    % % Plot log_10 sqrt (\int_x (U - U_true)^2 dx)
%     f3 = figure(13);
%     hold on;
%     L2_diff_U_true = sqrt(trapz(-L/2 : h : L/2, (U - true_sol).^2, 2));
%     pp = plot(linspace(0, T, length(L2_diff_U_true)), log10(L2_diff_U_true));
%     xl = xlabel('$t$ (s)', 'Interpreter', 'latex', 'FontSize', 18);
%     yl = ylabel('$\log_{10} \|u_j(t) - u_j^*(t)\|_2$', 'Interpreter', 'latex', 'FontSize', 18);
%     axis([0 T -5.5 -1]); pp.LineWidth = 1.4;
%     ax3 = f3.CurrentAxes; ax3.FontSize = 14;

% % Lazy animation
figure(14)
pause(1);
for i = 1 : 1 : Nt
    plot(0:h:L, U(i,:));
    axis([0 L h_min h_max]);
    title(['t = ' num2str((i-1)*tau)]);
    pause(0.01);
end
    
    
% %    Making animation
%     filename = strcat('abs_error_anim_', patern_left_str, 'nu_', num2str(nu), 'mu_', num2str(mu), '_tau_', num2str(tau), '_h_', num2str(h), '.gif');
%     figure(5);
%     pause(.5);
%     for i = 1 : ceil(.0005 / tau) : Nt
% %         plot(0:h:L, U(i, :));
% 
%         subplot(1,2,1);
%         plot(-L/2 : h : L/2, U(i, :)); hold on;
%         plot(-L/2 : h : L/2, true_sol(i, :), '--k'); hold off;
%         axis([-L/2 L/2 h_min h_max]);
%         xlabel('x', 'FontSize', 16);
%         lg = legend('u(t, x)', 'u_{true}(t, x)'); lg.FontSize = 14; lg.Location = 'northeast';
%         title({['t = ' num2str(tau * (i-1))], patern_left_str});
%         
%         subplot(1,2,2);
%         plot(-L/2 : h : L/2, log10(abs(U(i, :) - true_sol(i,:))));
%         title({['t = ' num2str(tau * (i-1))], patern_left_str});
%         axis([-L/2 L/2 -6 0]);
% %         axis([0 L -0.12 0.12]);
% %         ttl=title({['Волна в t = ', num2str((i - 1) * tau), 'с.'],...
% %             ['Шаблон: ', patern_str]}); ttl.FontSize = 16;
% %         ttl = title({'Волна в ', ['t = ', num2str((i - 1) * tau), 'с.']}); ttl.FontSize = 16;
%         xlabel('x', 'FontSize', 16);
%         yl = ylabel('log_{10} |\it{u - u^*}|'); yl.FontSize = 16;
% 
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