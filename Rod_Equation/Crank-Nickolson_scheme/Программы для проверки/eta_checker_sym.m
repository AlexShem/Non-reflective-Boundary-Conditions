nu = .1; mu = .15;
% nu = 1.4; mu = .5;
% nu = .01; mu = nu + 1/12 + .05;
a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*nu; s = nu/2;

syms w;

eq_eta = @(eta) s*(1+w.^2).*eta.^2 + (b*(1+w.^2)+g*w).*eta + a*(1+w.^2) + d*w - 2*s*(1+w.^2);

eps_1 = 1/(b^2-4*s*(a-2*s)) * (-b*g+2*d*s + 2*sqrt(-b*g*d*s+d^2*s^2+a*g^2*s-2*g^2*s^2));
eps_2 = 1/(b^2-4*s*(a-2*s)) * (-b*g+2*d*s - 2*sqrt(-b*g*d*s+d^2*s^2+a*g^2*s-2*g^2*s^2));

w_1 = .5*(eps_1 - sqrt(eps_1^2 - 4));
w_2 = .5*(eps_1 + sqrt(eps_1^2 - 4));
w_3 = .5*(eps_2 - sqrt(eps_2^2 - 4));
w_4 = .5*(eps_2 + sqrt(eps_2^2 - 4));

geom_sum = 0;
sum_1 = 0;
sum_2 = 0;
N = 101;
for n = 0 : N
    geom_sum = geom_sum + (-1)^n * w.^(2*n);
    sum_1 = sum_1 + P(n, eps_1/2) * w.^n;
    sum_2 = sum_2 + P(n, eps_2/2) * w.^n;
end

eta_1 = 1/(2*s) * geom_sum .* (-b*(1+w.^2) - g*w - sqrt(b^2 - 4*s*(a-2*s)) * (w.^2 - eps_1*w + 1) ...
    .* (w.^2 - eps_2*w + 1) .* sum_1 .* sum_2);
eta_2 = 1/(2*s) * geom_sum .* (-b*(1+w.^2) - g*w + sqrt(b^2 - 4*s*(a-2*s)) * (w.^2 - eps_1*w + 1) ...
    .* (w.^2 - eps_2*w + 1) .* sum_1 .* sum_2);

coef_1 = coeffs(expand(eta_1));
coef_1 = double(coef_1(1 : N));
coef_2 = coeffs(expand(eta_2));
coef_2 = double(coef_2(1 : N));

figure(2);
plot((1:N)-1, abs(coef_1)); hold on;
plot((1:N)-1, abs(coef_2), '.r');
xlabel('n');
legend('\eta_1', '\eta_2');
title({'Absolute value of \eta coefficients before \omega^n', ...
    ['\nu = ' num2str(nu) ', \mu = ' num2str(mu)]});