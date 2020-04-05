clear;
nu = .1;
mu = .2;
x = linspace(-1.2 , 1.2, 101); y = linspace(-1.2 , 1.2, 101);
[x, y] = meshgrid(x, y);
w = x + 1i*y;

a = 1 + 2*mu + 3*nu;
b = -mu - 2*nu;
g = 2*mu;
d = -2 - 4*mu;
s = nu/2;

eta_eq = @(eta) s*(1 + w.^2) .* eta.^2 + (b*(1 + w.^2) + g*w).*eta + d*w + (a- 2*s).*(1 + w.^2);

N = 350;
[eta_ser_minus, eta_ser_plus] = eta_series_num(nu, mu, N);

eta_minus = polyval(eta_ser_minus, w);
eta_plus = polyval(eta_ser_plus, w);

eq_minus = eta_eq(eta_minus);
eq_plus = eta_eq(eta_plus);

figure(1);
subplot(1,2,1);
contour(x, y, log10(abs(eq_minus)), 'ShowText', 'on');
xlabel('Real(\omega)'); ylabel('Imag(\omega)');
title('Equation for \eta_1(\omega)');
subplot(1,2,2);
contour(x, y, log10(abs(eq_plus)), 'ShowText', 'on');
xlabel('Real(\omega)'); ylabel('Imag(\omega)');
title('Equation for \eta_2(\omega)');