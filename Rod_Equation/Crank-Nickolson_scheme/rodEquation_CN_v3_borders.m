function U = rodEquation_CN_v3_borders(patterns)
%profile on;
% % Best [4 5 8 8], const condition; h = 0.02; tau = 1.6e-4; Steel
% parameters

%     h = 0.005; tau = 0.00001;
%     h = 0.005*2; tau = 0.00001*4;
h = 0.005*4; tau = 0.00001*16;

%     h = 0.005*4; tau = 0.00001*4;
%     h = 0.005*4; tau = 0.00001*1;

%     Set the parameters
%   par = [rho, R, E]
par = [7860, 1e-3, 210e9];  % %For Steel
%     par = [10000, 1e-2, 210e7];  % %RANDOM
%     par = [7860, 1.1e-3, 210e9];  % %Slight change
D = par(2)^2;   C = par(3) * par(2)^2 / par(1);

nu = C * tau^2 / h^4; % Параметр Куранта
mu = D / h^2;
a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*mu; s = nu/2;

L = 1;  %Длина струны, м.
%     T = .001; %Время, сек.
T = .3; %Время, сек.
%     T = .01*10; %Время, сек.

Nx = length(-L/2 : h : L/2); % Кол-во точек по пространству
Nt = length(0 : tau : T);   % Кол-во точек по времени

u_0 = linspace(0, L, Nx)';
%     u_0 = cosDistribFun_at_tao(0, 1, u_0);
x = linspace(-L/2, L/2, Nx); [u_0, sigma] = GaussianDistrib(-.3, .3, x); u_0 = u_0 .* x; u_0 = u_0(:);
h_min = min(u_0);
h_max = max(u_0);

    u_tau = u_0;
% %     f = GaussianDistrib(-.3, .3, (-20*L) : h : (20 * L)); f = f .* ((-20*L) : h : (20 * L)); f = f(:);
% f = GaussianDistrib(-.3, .3, linspace(-20*L, 20 * L, Nx)); f = f .* linspace(-20*L, 20 * L, Nx); f = f(:);
% u_tau = u_tau_function(h, D, C, f); len_tau = length(u_tau);
% u_tau = u_tau(ceil((len_tau+1)/2 - length(u_0)/2) : floor((len_tau+1)/2 + length(u_0)/2));
% u_tau = u_0 + u_tau * tau^2 / 2; u_tau = u_tau(:);

U = zeros(Nt, Nx);
U(1, :) = u_0';
U(2, :) = u_tau';
%     ---------------------------------------------------
%     polynomDEG_bor = [5, 5, 18, 18];
%     polynomDEG_prebor = [5, 5, 18, 18];

%     polynomDEG_bor = [5, 5, 10, 10];
%     polynomDEG_prebor = [5, 5, 10, 10];

%     polynomDEG_bor = [4, 4, 8, 8];
%     polynomDEG_prebor = [4, 4, 8, 8];

%     polynomDEG_bor = [3, 3, 9, 9];
%     polynomDEG_prebor = [3, 3, 9, 9];

%     polynomDEG_bor = [6, 4, 12, 8];
%     polynomDEG_prebor = [6, 4, 12, 8];

%     polynomDEG_bor = [5, 6, 5, 6];
%     polynomDEG_bor = [4, 5, 8, 8];

%    %Philipp work
%     polynomDEG_bor = [4, 5, 8, 8];
%     polynomDEG_bor = [5, 5, 10, 10];
%     polynomDEG_bor = [3, 5, 8, 8];
%     polynomDEG_bor_step = [5, 3, 8, 8];
%     polynomDEG_bor_step = [4, 4, 8, 8];
%     polynomDEG_bor = [5, 5, 11, 13];
%     polynomDEG_bor = [5, 7, 13, 16];
%     polynomDEG_bor = [6, 7, 13, 21];
%     polynomDEG_bor = [4, 4, 12, 11];
%     polynomDEG_bor = [5, 3, 3, 2];

%   % Seems good
polynomDEG_bor = [4, 5, 8, 8];
%     polynomDEG_bor = [6, 9, 7, 7];

polynomDEG_prebor = polynomDEG_bor;
%     ---------------------------------------------------
step_space_bor = length(polynomDEG_bor);
step_space_prebor = length(polynomDEG_prebor);

constant_condition = 1;
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

%     % Pure zero
    U_now([1, 2, end-1, end], [1, 2, end-1, end]) = eye(4);

%     % Zero and 0-second derivative
%     U_now(1, 1) = 1;
%     U_now(end, end) = 1;
%     U_now(2, 1) = 1; U_now(2, 2) = -2; U_now(2, 3) = 1;
%     U_now(end-1, end) = 0; U_now(end-1, end-1) = -2; U_now(end-1, end-2) = 1;

%     % Free rod
% U_now(1, 1) = 1; U_now(1, 3) = -3; U_now(1, 4) = 2;
% U_now(end, end) = 1; U_now(end, end-2) = -3; U_now(end, end-3) = 2;
% U_now(2, 2) = 1; U_now(2, 3) = -2; U_now(2, 4) = 1;
% U_now(end-1, end-1) = 1; U_now(end-1, end-2) = -2; U_now(end-1, end-3) = 1;

U_at_pTime = sparse(U_at_pTime);
U_at_ppTime = sparse(U_at_ppTime);
U_now = sparse(U_now);

%     % Main Cycle
for n = 3 : Nt
    rightPart = -U_at_pTime*U(n-1, :).'...
        -U_at_ppTime*U(n-2, :).';
    U(n, :) = U_now \ rightPart;
end

%profile report

%     patern_left = '';
%     patern_left_str = '';
%     for i = 1:length(polynomDEG_bor)
%         patern_left = strcat(patern_left, num2str(polynomDEG_bor(i)), '_');
%         patern_left_str = strcat(patern_left_str, num2str(polynomDEG_bor(i)), ',');
%     end
%     patern_left_str = patern_left_str(1:end-1);
%     patern_right = '';
%     patern_right_str = '';
%     for i = 1:length(polynomDEG_prebor)
%         patern_right = strcat(patern_right, num2str(polynomDEG_prebor(i)), '_');
%         patern_right_str = strcat(patern_right_str, num2str(polynomDEG_prebor(i)), ',');
%     end
%     patern_right_str = patern_right_str(1:end-1);

%     % True Solution
x_domain = linspace(0, L, Nx);
time = linspace(0, T, Nt);

true_sol = rodEquation_CN_PseudoTrueSolution_v2([h, tau], par, [nu mu], T, u_0, u_tau, x_domain);
center = size(true_sol, 2) / 2 + .5;
LB = ceil(center - Nx/2);
RB = floor(center + Nx/2);
true_sol = true_sol(:, LB:RB);

% % Plot log_10 |WaveEnergy - WaveEnergy_true|
%     figure(1); hold on;
%     I_minus_true = int_sum_abs2(U - true_sol, h);
%     pp = plot(linspace(0, T, length(I_minus_true)-1), log10(abs(I_minus_true(1:end-1)))); xl = xlabel('t, сек.');
% %     yl = ylabel('log_{10}I_d(t)');
% %     ttl = title({'log_{10}I_d(t)', ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)],...
% %         ['Bor: ', patern_str, ', PreBor: ', prepatern_str]});
% %     ttl = title({'log_{10} I(t)', ...
% %         ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)], ...
% %         [patern_left_str, ' ; ', patern_right_str]});
%     ttl = title({'log_{10} I(t)', ...
%         ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)], ...
%         ['\rho = ', num2str(par(1)), ', R =  ', num2str(par(2)), ', E = ', num2str(par(3))]});
%     axis([0 T -5.5 -1]); xl.FontSize = 18; yl.FontSize = 18; ttl.FontSize = 16; pp.LineWidth = 1.4;
% % %     ln = line([0.25 0.25], [-9.5 -2]); ln.LineStyle = '--'; ln.Color = 'k'; ln.LineWidth = 1.2;
%     grid on;

% % Plot Energy norm
    f1 = figure(11);
    hold on;
    H = int_sum_derivative(U - true_sol, h, tau, par(1), par(2), par(3));
    pp = plot(linspace(0, T, length(H)), log10(H));
    xl = xlabel('$t$ (s)', 'Interpreter', 'latex', 'FontSize', 18);
    yl = ylabel('$\log_{10} \mathcal{H} [u_j(t) - u_j^*(t)]$', 'Interpreter', 'latex', 'FontSize', 18);
    axis([0 T -1 4.5]); pp.LineWidth = 1.4;
    ax1 = f1.CurrentAxes; ax1.FontSize = 14;

% % Plot max log_10 |U - U_true|
    f2 = figure(12);
    hold on;
    Max_diff_U_true = max(abs(U - true_sol), [], 2);
    pp = plot(linspace(0, T, length(Max_diff_U_true)), log10(Max_diff_U_true));
    xl = xlabel('$t$ (s)', 'Interpreter', 'latex', 'FontSize', 18);
    yl = ylabel('$\log_{10} \max_j |u_j(t) - u_j^*(t)|$', 'Interpreter', 'latex', 'FontSize', 18);
    axis([0 T -5.5 -1]); pp.LineWidth = 1.4;
    ax2 = f2.CurrentAxes; ax2.FontSize = 14;
    
% % Plot log_10 sqrt (\int_x (U - U_true)^2 dx)
    f3 = figure(13);
    hold on;
    L2_diff_U_true = sqrt(trapz(-L/2 : h : L/2, (U - true_sol).^2, 2));
    pp = plot(linspace(0, T, length(L2_diff_U_true)), log10(L2_diff_U_true));
    xl = xlabel('$t$ (s)', 'Interpreter', 'latex', 'FontSize', 18);
    yl = ylabel('$\log_{10} \|u_j(t) - u_j^*(t)\|_2$', 'Interpreter', 'latex', 'FontSize', 18);
    axis([0 T -5.5 -1]); pp.LineWidth = 1.4;
    ax3 = f3.CurrentAxes; ax3.FontSize = 14;

% Lazy animation
% figure(6)
% pause(.5);
% for i = 1 : 2 : Nt
%     subplot(1,2,1);
%     plot(0:h:L, U(i,:)); hold on;
%     plot(0:h:L, true_sol(i, :), '--k'); hold off;
%     axis([0 L h_min h_max]);
%     title(['t = ' num2str((i-1)*tau)]);
%
%     subplot(1,2,2);
%     plot(0:h:L, log10(abs(U(i,:) - true_sol(i,:))));
%     title(['t = ' num2str((i-1)*tau)]);
%     axis([0 L -6 0])
%     pause(.01);
% end

end