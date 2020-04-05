nu = .1; mu = .01;

a = 3/(12*mu + 4) - .5; b = 9*nu/(3*mu + 1) - 2; c = 1 - (12*nu + 3) / (6*mu + 2); d = 3*nu / (6*mu + 2);

% w = linspace(1e-5, .05, 100);
% [x, y] = meshgrid(linspace(1e-5, .05, 60), linspace(-.05, .05, 60)); w = complex(x, y);
[x, y] = meshgrid(linspace(-.05, .05, 60), linspace(-.05, .05, 60)); w = complex(x, y);

% ETA Calculation
% eta_train = eta_true(nu, mu, w);
eta_plus = -(a*(1+w.^2) + c*w) ./ (2*d) + .5 * sqrt(((a*(1+w.^2)+c*w) ./ d).^2 - 4 * (w.^2 + (b-2*d)*w+1).*w ./ d);
eta_minus = -(a*(1+w.^2) + c*w) ./ (2*d) - .5 * sqrt(((a*(1+w.^2)+c*w) ./ d).^2 - 4 * (w.^2 + (b-2*d)*w+1).*w ./ d);

% GOOD
% eta_ser_plus = -(a*(1+w.^2)+c*w) ./ (2*d) + .5 * mysqrt((a^2*w.^4 + (2*a*c-4*d)*w.^3 ...
%     +(2*a^2+c^2-4*(b*d-2*d^2))*w.^2 + (2*a*c-4*d)*w + a^2) ./ d^2);

% GOOD
% eta_ser_plus = -(a*(1+w.^2)+c*w) ./ (2*d) + 1./(2*d) .* sqrt((a^2*w.^4 + (2*a*c-4*d)*w.^3 ...
%     +(2*a^2+c^2-4*(b*d-2*d^2))*w.^2 + (2*a*c-4*d)*w + a^2));

eps_2 = (2*d-a*c)/a^2 + 2/a^2 * sqrt(-a*c*d + d^2 - 2*a^2*d^2 + a^2*b*d);
w_1 = .5*(eps_2 - sqrt(eps_2^2-4));
w_2 = .5*(eps_2 + sqrt(eps_2^2-4));

% Quite good
% eta_ser_plus = -(a*(1+w.^2)+c*w) ./ (2*d) + 1./(2*d) .* sqrt(a^2 * (1-w).^2 .* (w-w_3).*(w-w_4));

% Quite good
% eta_ser_plus = -(a*(1+w.^2)+c*w) ./ (2*d) + a*(1-w)./(2*d) .* sqrt((w-w_1).*(w-w_2));

% eta_ser_plus = 1/(2*d) * (-a*(1+w.^2)-c*w + a*(1-w).* sqrt(w.^2 - (w_1 + w_2)*w + w_1*w_2));

% eta_ser_plus = 1/(2*d) * (-a*(1+w.^2)-c*w + a*(1-w) .* sqrt(w_1*w_2) .* sqrt(w.^2/(w_1*w_2) - (w_1 + w_2)/(w_1*w_2)*w + 1));

% eta_ser_plus = 1/(2*d) * (-a*(1+w.^2)-c*w + a*(1-w) .* sqrt(w.^2 - (w_1 + w_2)*w + 1));

N = 100;
eta_ser_plus = 0;
for n = 0 : N
    eta_ser_plus = eta_ser_plus + P(n, .5*(w_1+w_2)) * w.^n;
end
% eta_ser_minus = 1/(2*d) * (-a*(1+w.^2) - c*w - a*(1-w).*(w.^2-(w_1+w_2)*w+1).*eta_ser_plus); 
% eta_ser_plus = 1/(2*d) * (-a*(1+w.^2) - c*w + a*(1-w).*(w.^2-(w_1+w_2)*w+1).*eta_ser_plus);
eta_ser_minus = 1./(2*d*w) .* (-a*(1+w.^2) - c*w - a*(1-w).*(w.^2-(w_1+w_2)*w+1).*eta_ser_plus); 
eta_ser_plus = 1./(2*d*w) .* (-a*(1+w.^2) - c*w + a*(1-w).*(w.^2-(w_1+w_2)*w+1).*eta_ser_plus);

eq_minus = eta_ser_minus.^2 + (a*(1+w.^2)+c*w)./(d*w).*eta_ser_minus + (w.^2 + (b-2*d)*w + 1)./(d*w);
eq_plus = eta_ser_plus.^2 + (a*(1+w.^2)+c*w)./(d*w).*eta_ser_plus + (w.^2 + (b-2*d)*w + 1)./(d*w);

% figure(1);
% plot(log10(abs(w)), log10(abs(eta_plus - eta_ser_plus)));
% xlabel('\omega');
% title('log_{10} |\eta_+(\omega) - \eta_+^{series}(\omega)|');
% 
% figure(2);
% plot(log10(abs(w)), log10(abs(eta_minus - eta_ser_minus)));
% xlabel('\omega');
% title('log_{10} |\eta_-(\omega) - \eta_-^{series}(\omega)|');

% figure(3);
% contour(real(w), imag(w), log10(abs(eta_ser_plus - eta_plus)), 'ShowText', 'on');
% xlabel('Re(\omega)'); ylabel('Im(\omega)');
% title('log_{10} |\eta_1(\omega) - \eta_1^{series}(\omega)|');
% 
% figure(4);
% contour(real(w), imag(w), log10(abs(eta_ser_minus - eta_minus)), 'ShowText', 'on');
% xlabel('Re(\omega)'); ylabel('Im(\omega)');
% title('log_{10} |\eta_2(\omega) - \eta_2^{series}(\omega)|');

figure(3);
contour(real(w), imag(w), log10(abs(eq_minus)), 'ShowText', 'on');
xlabel('Re(\omega)'); ylabel('Im(\omega)');
title({'log_{10} |\eta_-^2(\omega) - ... \eta_-(\omega) + ...|',...
    ['N_{Leg} = ' num2str(N)]});

figure(4);
contour(real(w), imag(w), log10(abs(eq_plus)), 'ShowText', 'on');
xlabel('Re(\omega)'); ylabel('Im(\omega)');
title({'log_{10} |\eta_+^2(\omega) - ... \eta_+(\omega) + ...|',...
    ['N_{Leg} = ' num2str(N)]});