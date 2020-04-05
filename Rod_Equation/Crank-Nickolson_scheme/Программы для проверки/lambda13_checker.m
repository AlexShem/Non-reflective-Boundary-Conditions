nu = .1; mu = .15;
% nu = .2441; mu = .0001;
% nu = .01; mu = nu + 1/12 + .05;
a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*nu; s = nu/2;

[x, y] = meshgrid(linspace(-.7, .7, 81), linspace(-.7, .7, 81));
% [x, y] = meshgrid(linspace(-.1, .1, 81), linspace(-.1, .1, 81));
w = x+1i*y;

% w = linspace(1e-8, .01, 1001);

N_eta = 100;
[eta_1, ~] = eta_fun(nu, mu, w, N_eta);

eq_lam = @(lam) lam.^2 - eta_1 .* lam + 1;

theta = 1/(2*s) * (-b - sqrt(b^2 - 4*s*(a-2*s)));
r = eta_1 - theta;

% lam_1 = eta_1/2 - sqrt(.25 * (theta + r).^2 - 1);
% lam_3 = eta_1/2 + sqrt(.25 * (theta + r).^2 - 1);

% lam_1 = eta_1/2 - sqrt(.25 * (r.^2 + 2*theta*r + theta^2) - 1);
% lam_3 = eta_1/2 + sqrt(.25 * (r.^2 + 2*theta*r + theta^2) - 1);

% lam_1 = eta_1/2 - sqrt(theta^2/4 - 1) * sqrt(r.^2 / (theta^2 - 4) + 2*theta*r / (theta^2 - 4) + 1);
% lam_3 = eta_1/2 + sqrt(theta^2/4 - 1) * sqrt(r.^2 / (theta^2 - 4) + 2*theta*r / (theta^2 - 4) + 1);

zeta = 1/sqrt(theta^2 - 4) * r;

% lam_1 = eta_1/2 - sqrt(theta^2/4 - 1) * sqrt(zeta.^2 + 2*theta*zeta / sqrt(theta^2 - 4) + 1);
% lam_3 = eta_1/2 + sqrt(theta^2/4 - 1) * sqrt(zeta.^2 + 2*theta*zeta / sqrt(theta^2 - 4) + 1);

lam_1 = eta_1/2 - sqrt(theta^2/4 - 1) * (zeta.^2 + 2*theta*zeta / sqrt(theta^2 - 4) + 1) .* (zeta.^2 + 2*theta*zeta / sqrt(theta^2 - 4) + 1).^-.5;
lam_3 = eta_1/2 + sqrt(theta^2/4 - 1) * (zeta.^2 + 2*theta*zeta / sqrt(theta^2 - 4) + 1) .* (zeta.^2 + 2*theta*zeta / sqrt(theta^2 - 4) + 1).^-.5;

% N_lam = N_eta;
% sum = 0;
% for n = 0 : N_lam
%     sum = sum + P(n, -theta/sqrt(theta^2-4)) * zeta.^n;
% end
% lam_1 = eta_1/2 - sqrt(theta^2/4 - 1) * (zeta.^2 + 2*theta*zeta / sqrt(theta^2 - 4) + 1) .* sum;
% lam_3 = eta_1/2 + sqrt(theta^2/4 - 1) * (zeta.^2 + 2*theta*zeta / sqrt(theta^2 - 4) + 1) .* sum;

% N_lam = N_eta;
% sum = 0;
% for n = 0 : N_lam
%     sum = sum + P(n, -theta/sqrt(theta^2-4)) * (theta^2 - 4)^(n/2) * r.^n;
% end
% lam_1 = eta_1/2 - sqrt(theta^2/4 - 1) * (1/(theta^2-4)*r.^2 + 2*theta/(theta^2-4)*r / (theta^2 - 4) + 1) .* sum;
% lam_3 = eta_1/2 + sqrt(theta^2/4 - 1) * (1/(theta^2-4)*r.^2 + 2*theta/(theta^2-4)*r / (theta^2 - 4) + 1) .* sum;

eq_1 = eq_lam(lam_1);
eq_3 = eq_lam(lam_3);

figure(1);
contour(x,y, log10(abs(eq_1)), 'ShowText', 'on');
title('log_{10} |\lambda_1^2 - \eta_1\lambda_1 + 1|');
xlabel('\nu'); ylabel('\mu');

figure(2);
contour(x,y, log10(abs(eq_3)), 'ShowText', 'on');
title('log_{10} |\lambda_3^2 - \eta_1\lambda_3 + 1|');
xlabel('\nu'); ylabel('\mu');