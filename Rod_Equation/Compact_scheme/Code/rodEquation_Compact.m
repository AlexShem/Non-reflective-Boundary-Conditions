function U = rodEquation_Compact(patterns)
%profile on

%     h = 0.005; tau = 0.001;   %Precise
h = 0.0025*8; tau = 0.00001;
%     h = 0.00025; tau = 0.00025;
%     h = 0.001; tau = 0.0002;

%     Set the parameters
%   par = [rho, R, E]
par = [7860, 1e-2, 210e9];  % %For Steel
D = par(2)^2;   C = par(3) * par(2)^2 / par(1);

nu = C * tau^2 / h^4;
mu = D / h^2;
a = 3/(12*mu + 4) - .5; b = 9*nu/(3*mu + 1) - 2; c = 1 - (12*nu + 3) / (6*mu + 2); d = 3*nu / (6*mu + 2);

L = 1;
T = 1;

Nx = single(L / h + 1);
Nt = single(T / tau + 1);

x = linspace(0, L, Nx);
time = linspace(0, T, Nt);
[x, time] = meshgrid(x, time);
true_sol = wave_true_sol(x, time, D) * .1;

u_0 = linspace(0, L, Nx)';
x = -L/2 : h : L/2; [u_0, sigma] = GaussianDistrib(-.3, .3, x); u_0 = u_0 .* x; u_0 = u_0(:);
h_min = min(u_0);
h_max = max(u_0);

%     u_tau = u_0;
f = GaussianDistrib(-.3, .3, (-20*L) : h : (20 * L)); f = f .* ((-20*L) : h : (20 * L)); f = f(:);
u_tau = u_tau_function(h, D, C, f); len_tau = length(u_tau);
u_tau = u_tau(ceil((len_tau+1)/2 - length(u_0)/2) : floor((len_tau+1)/2 + length(u_0)/2));
u_tau = u_0 + u_tau * tau^2 / 2; u_tau = u_tau(:);

U = zeros(Nt, Nx);
U(1, :) = u_0';
U(2, :) = u_tau';
%     ---------------------------------------------------
%     polynomDEG_bor = [0, 0, 2, 2];
%     polynomDEG_prebor = [0, 0, 2, 2];
%     polynomDEG_bor = [0, 0, 3, 3];
%     polynomDEG_prebor = [0, 0, 3, 3];

polynomDEG_bor = [4, 5, 8, 8];


polynomDEG_prebor = polynomDEG_bor;
%     ---------------------------------------------------
step_space_bor = length(polynomDEG_bor);
step_space_prebor = length(polynomDEG_prebor);

constant_condition = 0;
x_condition = 0;

[borderDEG, preborderDEG] = ...
    PadeCoef_poly_const_x(nu, mu, polynomDEG_bor, polynomDEG_prebor, ...
    constant_condition, x_condition);

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

U_now = diag(a*ones(1, Nx)) + diag(b * ones(1, Nx-1), 1) + diag(b * ones(1, Nx-1), -1) + diag(s * ones(1, Nx-2), 2) + diag(s * ones(1, Nx-2), -2);
U_now(1, :) = 0; U_now(2, :) = 0;
U_now(end, :) = 0; U_now(end-1, :) = 0;

U_at_pTime = diag(d * ones(1, Nx)) + diag(g * ones(1, Nx-1), 1) + diag(g * ones(1, Nx-1), -1);
U_at_pTime(1, :) = 0; U_at_pTime(2, :) = 0;
U_at_pTime(end, :) = 0; U_at_pTime(end-1, :) = 0;

U_at_ppTime = diag(a*ones(1, Nx)) + diag(b * ones(1, Nx-1), 1) + diag(b * ones(1, Nx-1), -1) + diag(s * ones(1, Nx-2), 2) + diag(s * ones(1, Nx-2), -2);
U_at_ppTime(1, :) = 0; U_at_ppTime(2, :) = 0;
U_at_ppTime(end, :) = 0; U_at_ppTime(end-1, :) = 0;

U_now(1:2, 1:step_space_bor) = [coefsm_bor(1, :); coefsm_prebor(1, :)];
U_now(end:-1:end-1, end:-1:end-step_space_bor+1) = [coefsm_bor(1, :); coefsm_prebor(1, :)];
U_at_pTime(1:2, 1:step_space_bor) = [coefsm_bor(2, :); coefsm_prebor(2, :)];
U_at_pTime(end:-1:end-1, end:-1:end-step_space_bor+1) = [coefsm_bor(2, :); coefsm_prebor(2, :)];
U_at_ppTime(1:2, 1:step_space_bor) = [coefsm_bor(3, :); coefsm_prebor(3, :)];
U_at_ppTime(end:-1:end-1, end:-1:end-step_space_bor+1) = [coefsm_bor(3, :); coefsm_prebor(3, :)];



U_at_pTime = sparse(U_at_pTime);
U_at_ppTime = sparse(U_at_ppTime);
U_now = sparse(U_now);
if max(abs(condest(U_now))) > 1e9
    warning('Matrix U_now is degenerate');
end
U_now_inv = inv(U_now);

nStepT_bor = max(polynomDEG_bor) + 1;
nStepT_prebor = max(polynomDEG_prebor) + 1;

for n = 3 : Nt
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
    
    U(n, :) = - U_now_inv * rightPart;
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

%     % True Solution
x_domain = linspace(0, L, Nx);
time = linspace(0, T, Nt);

true_sol = rodEquation_CN_PseudoTrueSolution([h, tau], [nu mu], T, u_0, u_tau, x_domain);
center = size(true_sol, 2) / 2 + .5;
LB = ceil(center - Nx/2);
RB = floor(center + Nx/2);
true_sol = true_sol(:, LB:RB);


% % Plot max log_10 |U - U_true|
f2 = figure(2);
hold on;
Max_diff_U_true = max(abs(U - true_sol), [], 2);
pp = plot(linspace(0, T, length(Max_diff_U_true)), log10(Max_diff_U_true)); xl = xlabel('t');
%     yl = ylabel('log_{10}M_d(t)');
%     ttl = title({'log_{10}M_d(t)', ['\nu = ', num2str(nu) '; h = ', num2str(h), '; c = ', num2str(c)], ['Шаблон: ', patern_str]});
%     ttl = title({'log_{10} max_x |\hat{u} - u| (t)', ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)],...
%         ['Bor: ', patern_str, ', PreBor: ', prepatern_str]});
%     ttl = title({'log_{10} |u - u^*|', ...
%         ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)], ...
%         [patern_left_str, ' ; ', patern_right_str]});
%     ttl = title({'log_{10} |u - u^*|', ...
%         ['\nu = ', num2str(nu), ', \mu = ', num2str(mu), ', h = ', num2str(h), ', \tau = ', num2str(tau)], ...
%         ['\rho = ', num2str(par(1)), ', R =  ', num2str(par(2)), ', E = ', num2str(par(3))]});
%     ttl.FontSize = 16;
axis([0 T -5.5 -1]); xl.FontSize = 18; pp.LineWidth = 1.4;
% %     ln = line([0.25 0.25], [-6 -2.5]); ln.LineStyle = '--'; ln.Color = 'k'; ln.LineWidth = 1.2;
grid on;
ax2 = f2.CurrentAxes; ax2.FontSize = 14;
yl = ylabel('max_j lg |u_j - u_j^*|');
yl.FontSize = 18;



% Plot max log_10 |U - U_true|
f3 = figure(3);
hold on;
relative_error = max(log10(abs(U - true_sol)), [], 2) ./ log10(max(true_sol, [], 2));
%     relative_error = max((abs(U - true_sol)), [], 2) ./ (max(true_sol, [], 2));
pp = plot(linspace(0, T, length(relative_error)), relative_error); xl = xlabel('t');
axis([0 T 0 8]); xl.FontSize = 18; pp.LineWidth = 1.4;
grid on;
ax3 = f3.CurrentAxes; ax3.FontSize = 14;
yl = ylabel('max_j lg |u_j - u_j^*| / lg max_j |u_j^*|');
yl.FontSize = 18;


% Lazy animation
% figure(6)
% pause(.1);
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


% %    Making animation
% filename = strcat('abs_error_anim_', patern_left_str, 'nu_', num2str(nu), 'mu_', num2str(mu), '_tau_', num2str(tau), '_h_', num2str(h), '.gif');
% figure(5);
% pause(.5);
% for i = 1 : ceil(.0005 / tau) : Nt
%     plot(0:h:L, U(i, :));
%
%     subplot(1,2,1);
%     plot(-L/2 : h : L/2, U(i, :)); hold on;
%     plot(-L/2 : h : L/2, true_sol(i, :), '--k'); hold off;
%     axis([-L/2 L/2 h_min h_max]);
%     xlabel('x', 'FontSize', 16);
%     lg = legend('u(t, x)', 'u_{true}(t, x)'); lg.FontSize = 14; lg.Location = 'northeast';
%     title({['t = ' num2str(tau * (i-1))], patern_left_str});
%
%     subplot(1,2,2);
%     plot(-L/2 : h : L/2, log10(abs(U(i, :) - true_sol(i,:))));
%     title({['t = ' num2str(tau * (i-1))], patern_left_str});
%     axis([-L/2 L/2 -6 0]);
%     axis([0 L -0.12 0.12]);
%     ttl=title({['Волна в t = ', num2str((i - 1) * tau), 'с.'],...
%         ['Шаблон: ', patern_str]}); ttl.FontSize = 16;
%     ttl = title({'Волна в ', ['t = ', num2str((i - 1) * tau), 'с.']}); ttl.FontSize = 16;
%     xlabel('x', 'FontSize', 16);
%     yl = ylabel('log_{10} |\it{u - u^*}|'); yl.FontSize = 16;
%
%     drawnow;
%     frame = getframe(5);
%     im = frame2im(frame);
%     [A, map] = rgb2ind(im, 256);
%     if i == 1
%         imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.001);
%     else
%         imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0);
%     end
% end
end