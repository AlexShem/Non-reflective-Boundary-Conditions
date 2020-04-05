function U = rodEquation_CN_htau_search(h, tau, polynomDEG_bor, T)
%     Set the parameters
par = [7860, 1e-3, 210e9];  % %For Steel
%     par = [7860, 1.1e-3, 210e9];  % %Slight change
D = par(2)^2;   C = par(3) * par(2)^2 / par(1);

nu = C * tau^2 / h^4;
mu = D / h^2;
a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*mu; s = nu/2;

L = 1;
% T = .04*8;

Nx = floor(L / h + 1);
Nt = floor(T / tau + 1);

u_0 = linspace(0, L, Nx)';
x = -L/2 : h : L/2;
[u_0, sigma] = GaussianDistrib(-.3, .3, x);
u_0 = u_0 .* x; u_0 = u_0(:);
h_min = min(u_0);
h_max = max(u_0);

%     u_tau = u_0;
% f = GaussianDistrib(-.3, .3, (-20*L) : h : (20*L)); f = f .* ((-20*L) : h : (20*L)); f = f(:);
% u_tau = u_tau_function(h, D, C, f); len_tau = length(u_tau);
% u_tau = u_tau(ceil((len_tau+1)/2 - length(u_0)/2) : floor((len_tau+1)/2 + length(u_0)/2));
% u_tau = u_0 + u_tau * tau^2 / 2; u_tau = u_tau(:);
u_tau = u_0;

U = zeros(Nt, Nx);
U(1, :) = u_0';
U(2, :) = u_tau';
%     ---------------------------------------------------
polynomDEG_prebor = polynomDEG_bor;
%     ---------------------------------------------------
step_space_bor = length(polynomDEG_bor);
step_space_prebor = length(polynomDEG_prebor);

if mod(sum(polynomDEG_bor), 2) == 1
    constant_condition = 1;
else
    constant_condition = 0;
end

[borderDEG, preborderDEG] = ...
    PadeCoef_poly_const_x(nu, mu, polynomDEG_bor, polynomDEG_prebor, ...
    constant_condition, 0);

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


U_now = spdiags(ones(Nx, 1) .* [s b a b s], -2:2, Nx, Nx);
U_at_pTime = spdiags(ones(Nx, 1) .* [g d g], -1:1, Nx, Nx);
U_at_ppTime = U_now;

U_now(1:2, 1:step_space_bor) = [coefsm_bor(1, :); coefsm_prebor(1, :)];
U_now(end:-1:end-1, end:-1:end-step_space_bor+1) = [coefsm_bor(1, :); coefsm_prebor(1, :)];
U_at_pTime(1:2, 1:step_space_bor) = [coefsm_bor(2, :); coefsm_prebor(2, :)];
U_at_pTime(end:-1:end-1, end:-1:end-step_space_bor+1) = [coefsm_bor(2, :); coefsm_prebor(2, :)];
U_at_ppTime(1:2, 1:step_space_bor) = [coefsm_bor(3, :); coefsm_prebor(3, :)];
U_at_ppTime(end:-1:end-1, end:-1:end-step_space_bor+1) = [coefsm_bor(3, :); coefsm_prebor(3, :)];

nStepT_bor = max(polynomDEG_bor) + 1;
nStepT_prebor = max(polynomDEG_prebor) + 1;

%     % Main Cycle
contr = max(abs(u_0));
for n = 3 : Nt
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
    
    if (max(abs(U(n, :)))) > 100*contr
        error('Scheme explodes');
    end
end

% x_domain = linspace(0, L, Nx);
% time = linspace(0, T, Nt);

% tic
% true_sol = rodEquation_CN_PseudoTrueSolution_v2([h, tau], par, [nu mu], T, u_0, u_tau, x_domain);
% toc
% center = size(true_sol, 2) / 2 + .5;
% LB = ceil(center - Nx/2);
% RB = floor(center + Nx/2);
% true_sol = true_sol(:, LB:RB);

% % Lazy animation
% figure(6)
% pause(1);
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