%% Intro
addpath('..\Code')

nu = 5;
mu = 1;

a = 3/(12*mu + 4) - .5;
b = (9*nu)/(3*mu + 1) - 2;
c = 1 - (12*nu + 3)/(6*mu + 2);
d = (3*nu)/(6*mu + 2);

x = linspace(-.1, .1, 1000);
y = linspace(-.1, .1, 1000);
[x, y] = meshgrid(x, y);
w = complex(x, y);



eta_min = (-a*(1+w.^2) - c*w - sqrt((a*(1+w.^2) + c*w).^2 - 4*d*w.*(w.^2 + (b-2*d)*w + 1))) ./ (2*d*w);
eta_pl = (-a*(1+w.^2) - c*w + sqrt((a*(1+w.^2) + c*w).^2 - 4*d*w.*(w.^2 + (b-2*d)*w + 1))) ./ (2*d*w);



eq_eta = @(eta) d*w*eta.^2 + (a*(1+w.^2) + c*w).*eta + w.^2 + (b-2*d)*w + 1;

% Chech eta_min
eq_eta_min_val = eq_eta(eta_min);
figure(1)
contour(x, y, log(abs(eq_eta_min_val)), 'ShowText', 'on');
% contour(x, y, abs(eta_min), 'ShowText', 'on');
title('Eq(\eta_{-}(\omega)) ?= 0');

% Chech eta_pl
eq_eta_pl_val = eq_eta(eta_pl);
figure(2)
contour(x, y, log(abs(eq_eta_pl_val)), 'ShowText', 'on');
% contour(x, y, abs(eta_pl), 'ShowText', 'on');
title('Eq(\eta_{+}(\omega)) ?= 0')



%% Series
N_terms = 20;
[eta_ser_minus, eta_ser_plus] = series_eta(nu, mu, N_terms);

eta_ser_minus_val = polyval(eta_ser_minus, w);
% Chech eta_min
eq_eta_min_val = eq_eta(eta_ser_minus_val);
figure(3)
contour(x, y, log(abs(eq_eta_min_val)), 'ShowText', 'on');
% contour(x, y, abs(eta_min), 'ShowText', 'on');
title('\eta_{-}(\omega)');