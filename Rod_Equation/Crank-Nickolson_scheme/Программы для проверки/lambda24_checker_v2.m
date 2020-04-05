nu = .1; mu = .15;
% nu = .2441; mu = .0001;
% nu = .01; mu = nu + 1/12 + .05;
a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*nu; s = nu/2;

[x, y] = meshgrid(linspace(-.7, .7, 81), linspace(-.7, .7, 81));
% [x, y] = meshgrid(linspace(-.1, .1, 81), linspace(-.1, .1, 81));
w = x+1i*y;

% w = linspace(1e-8, .01, 1001);

N_eta = 100;
[~, eta_2] = eta_fun(nu, mu, w, N_eta);

eq_lam = @(lam) lam.^2 - eta_2 .* lam + 1;
eq_lam_full = @(lam) s*(w.^2+1).*(lam.^4+1) + (b*(w.^2+1)+g*w).*(lam.^3+lam) + (a*(1+w.^2)+d*w).*lam.^2;

theta = 1/(2*s) * (-b + sqrt(b^2 - 4*s*(a-2*s)));
r = eta_2 - theta;

% lam_1 = eta_1/2 - sqrt(.25 * theta^2 - 1) * sqrt(1 + r/(theta + 2)) .* sqrt(1 + r/(theta - 2));
% lam_3 = eta_1/2 - sqrt(.25 * theta^2 - 1) * sqrt(1 + r/(theta + 2)) .* sqrt(1 + r/(theta - 2));

N_sum = 40;
sum_plus = 0;
sum_minus = 0;

mul = 1;
mul_plus = 1;
mul_minus = 1;

for n = 0 : N_sum
    sum_plus = sum_plus + mul * mul_plus * r.^n;
    sum_minus = sum_minus + mul * mul_minus * r.^n;
    
    mul = -mul * (2*n+1) * (2*n+2) * (1-2*n) / ((-1 - 2*n) * (n+1)^2 * 4);
    mul_plus = mul_plus / (theta + 2);
    mul_minus = mul_minus / (theta - 2);
end

% for n = 0 : N_sum
%     sum_plus = sum_plus + (-1)^n * factorial(2*n) / ((1-2*n) * factorial(n)^2 * 4^n * (theta + 2)^n) * r.^n;
%     sum_minus = sum_minus + (-1)^n * factorial(2*n) / ((1-2*n) * factorial(n)^2 * 4^n * (theta - 2)^n) * r.^n;
% end

lam_2 = eta_2/2 - sqrt(.25 * theta^2 - 1) * sum_plus .* sum_minus;
lam_4 = eta_2/2 + sqrt(.25 * theta^2 - 1) * sum_plus .* sum_minus;

% eq_2 = eq_lam(lam_2);
% eq_4 = eq_lam(lam_4);
eq_2 = eq_lam_full(lam_2);
eq_4 = eq_lam_full(lam_4);

figure(1);
contour(x,y, log10(abs(eq_2)), 'ShowText', 'on');
title({'log_{10} |\lambda_2^2 - \eta_2\lambda_2 + 1|', ['N = ' num2str(N_sum), ', \nu = ', num2str(nu), ', \mu = ', num2str(mu)]});
xlabel('Real(\omega)'); ylabel('Imag(\omega)');

figure(2);
contour(x,y, log10(abs(eq_4)), 'ShowText', 'on');
title({'log_{10} |\lambda_4^2 - \eta_2\lambda_4 + 1|', ['N = ' num2str(N_sum), ', \nu = ', num2str(nu), ', \mu = ', num2str(mu)]});
xlabel('Real(\omega)'); ylabel('Imag(\omega)');