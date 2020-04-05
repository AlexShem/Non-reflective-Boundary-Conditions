% nu = .1; mu = .01;
nu = 20100; mu = 10000;

a = 3/(12*mu + 4) - .5; b = 9*nu/(3*mu + 1) - 2; c = 1 - (12*nu + 3) / (6*mu + 2); d = 3*nu / (6*mu + 2);

% w = linspace(1e-5, .05, 100);
% [x, y] = meshgrid(linspace(1e-5, .05, 60), linspace(-.05, .05, 60)); w = complex(x, y);
[x, y] = meshgrid(linspace(-.05, .05, 60), linspace(-.05, .05, 60)); w = complex(x, y);

% eta_1 = -(a*(1+w.^2) + c*w) ./ (2*d*w) + .5 * sqrt(((a*(1+w.^2)+c*w) ./ (d*w)).^2 - 4 * (w.^2 + (b-2*d)*w+1) ./ (d*w));
% % eta_2 = -(a*(1+w.^2) + c*w) ./ (2*d*w) - .5 * sqrt(((a*(1+w.^2)+c*w) ./ (d*w)).^2 - 4 * (w.^2 + (b-2*d)*w+1) ./ (d*w));
% % eta_1 = -(a*(1+w.^2) + c*w) ./ (2*d) + .5 * sqrt(((a*(1+w.^2)+c*w) ./ (d)).^2 - 4 * (w.^2 + (b-2*d)*w+1).*w ./ (d));
% eta_2 = -(a*(1+w.^2) + c*w) ./ (2*d*w) - .5 * sqrt(((a*(1+w.^2)+c*w) ./ (d*w)).^2 - 4 * (w.^2 + (b-2*d)*w+1) ./ (d*w));

n_eta = 15;
[eta_1, eta_2] = eta_series(w, nu, mu, n_eta);

lam_1 = eta_1 / 2 - sqrt(eta_1.^2 / 4 - 1);
lam_2 = eta_2 / 2 - sqrt(eta_2.^2 / 4 - 1);
lam_3 = eta_1 / 2 + sqrt(eta_1.^2 / 4 - 1);
lam_4 = eta_2 / 2 + sqrt(eta_2.^2 / 4 - 1);

r = eta_1 + 1/a;
% lam_ser_1 = eta_1 / 2 - sqrt(1/(4*a^2)-1) * sqrt(a^2/(1-4*a^2)*r.^2 - 2*a/(1-4*a^2)*r + 1);
% lam_ser_3 = eta_1 / 2 + sqrt(1/(4*a^2)-1) * sqrt(a^2/(1-4*a^2)*r.^2 - 2*a/(1-4*a^2)*r + 1);

zeta = a/sqrt(1-4*a^2) * r;
% lam_ser_1 = eta_1/2 - sqrt(1/(4*a^2)-1) * sqrt(zeta.^2 - 2/sqrt(1-4*a^2)*zeta + 1);
% lam_ser_3 = eta_1/2 + sqrt(1/(4*a^2)-1) * sqrt(zeta.^2 - 2/sqrt(1-4*a^2)*zeta + 1);

% N = 3;
% lam_ser_1 = 0;
% for n = 0 : N
%     lam_ser_1 = lam_ser_1 + P(n, 1/sqrt(1-4*a^2)) * zeta.^n;
% end
% lam_ser_3 = eta_1/2 + sqrt(1/(4*a^2)-1) * (zeta.^2 - 2/sqrt(1-4*a^2)*zeta + 1) .* lam_ser_1;
% lam_ser_1 = eta_1/2 - sqrt(1/(4*a^2)-1) * (zeta.^2 - 2/sqrt(1-4*a^2)*zeta + 1) .* lam_ser_1;

% lam_2 and lam_4
w_1 = (-b + 2*d - sqrt((b-2*d)^2-4))/2;
w_2 = (-b + 2*d + sqrt((b-2*d)^2-4))/2;
A = -1/sqrt((b-2*d)^2-4);
fraction = 0;
                                N_frac = 10;
for k = 0 : N_frac
    fraction = fraction + (1/w_2^(k+1) - 1/w_1^(k+1)) * w.^k;
end
lam_ser_2 = 0;
                                N_lam2_4 = 10;
mul_fac = 1;
mul = 1;
for n = 0 : N_lam2_4
    lam_ser_2 = lam_ser_2 + mul * mul_fac / (1-2*n) * (w .* eta_1 .* fraction).^(2*n);
    mul = mul * (A*d)^2;
    mul_fac = mul_fac * (2*n+1)*(2*n+2)/(n+1)^2;
end
lam_ser_4 = eta_2/2 .* (1 + lam_ser_2);
lam_ser_2 = eta_2/2 .* (1 - lam_ser_2);
eq_2 = lam_ser_2.^2 - eta_2 .* lam_ser_2 + 1;
eq_4 = lam_ser_4.^2 - eta_2 .* lam_ser_4 + 1;

N = 50;
lam_ser_1 = 0;
for n = 0 : N
    lam_ser_1 = lam_ser_1 + P(n, 1/sqrt(1-4*a^2)) * (a/sqrt(1-4*a^2))^n * r.^n;
end
lam_ser_3 = eta_1/2 + sqrt(1/(4*a^2)-1) * (a^2/(1-4*a^2)*r.^2 - 2*a/(1-4*a^2)*r + 1) .* lam_ser_1;
lam_ser_1 = eta_1/2 - sqrt(1/(4*a^2)-1) * (a^2/(1-4*a^2)*r.^2 - 2*a/(1-4*a^2)*r + 1) .* lam_ser_1;

eq_1 = lam_ser_1.^2 - eta_1 .* lam_ser_1 + 1;
eq_3 = lam_ser_3.^2 - eta_1 .* lam_ser_3 + 1;

% figure(1);
% plot(log10(w), log10(abs(lam_ser_1 - lam_1)));
% xlabel('log_{10} \omega');
% title('log_{10} |\lambda_1(\omega) - \lambda_1^{series}(\omega)|');
% 
% figure(2);
% plot(log10(w), log10(abs(lam_ser_3 - lam_3)));
% xlabel('log_{10} \omega');
% title('log_{10} |\lambda_3(\omega) - \lambda_3^{series}(\omega)|');

% figure(3);
% contour(real(w), imag(w), log10(abs(lam_ser_1 - lam_1)), 'ShowText', 'on');
% xlabel('Re(\omega)'); ylabel('Im(\omega)');
% title('log_{10} |\lambda_1(\omega) - \lambda_1^{series}(\omega)|');
% 
% figure(4);
% contour(real(w), imag(w), log10(abs(lam_ser_3 - lam_3)), 'ShowText', 'on');
% xlabel('Re(\omega)'); ylabel('Im(\omega)');
% title('log_{10} |\lambda_3(\omega) - \lambda_3^{series}(\omega)|');

figure(3);
contour(real(w), imag(w), log10(abs(eq_1)), 'ShowText', 'on');
xlabel('Re(\omega)'); ylabel('Im(\omega)');
title({'log_{10} |\lambda_1^2(\omega) - \eta_1(\omega) \lambda_1(\omega) + 1|',...
    ['N_{Leg} = ' num2str(N) ', N_{\eta} = ' num2str(n_eta)]});

figure(4);
contour(real(w), imag(w), log10(abs(eq_3)), 'ShowText', 'on');
xlabel('Re(\omega)'); ylabel('Im(\omega)');
title({'log_{10} |\lambda_3^2(\omega) - \eta_1(\omega) \lambda_3(\omega) + 1|',...
    ['N_{Leg} = ' num2str(N) ', N_{\eta} = ' num2str(n_eta)]});

figure(5);
contour(real(w), imag(w), log10(abs(eq_2)), 'ShowText', 'on');
xlabel('Re(\omega)'); ylabel('Im(\omega)');
title({'log_{10} |\lambda_2^2(\omega) - \eta_2(\omega) \lambda_2(\omega) + 1|',...
    ['N_{Lam} = ' num2str(N_lam2_4) ', N_{Fraction} = ' num2str(N_frac) ', N_{\eta} = ' num2str(n_eta)]});

figure(6);
contour(real(w), imag(w), log10(abs(eq_4)), 'ShowText', 'on');
xlabel('Re(\omega)'); ylabel('Im(\omega)');
title({'log_{10} |\lambda_4^2(\omega) - \eta_2(\omega) \lambda_4(\omega) + 1|',...
    ['N_{Lam} = ' num2str(N_lam2_4) ', N_{Fraction} = ' num2str(N_frac) ', N_{\eta} = ' num2str(n_eta)]});