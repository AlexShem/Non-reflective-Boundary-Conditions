clear;
nu = .8;
mu = .1;
x = linspace(-1.2 , 1.2, 101); y = linspace(-1.2 , 1.2, 101);
[x, y] = meshgrid(x, y);
w = x + 1i*y;

a = 1 + 2*mu + 3*nu;
b = -mu - 2*nu;
g = 2*mu;
d = -2 - 4*mu;
s = nu/2;

lam_eq = @(lam) s*(1 + w.^2) .* (lam.^4 + 1) + (b*(1 + w.^2) + g*w) .* (lam.^3 + lam) + (a*(1 + w.^2) + d*w).*lam.^2;

N = 255;
% [lam_ser_1, lam_ser_2, lam_ser_3, lam_ser_4] = lambda_series_num(nu, mu, N);
[lam_ser_1, lam_ser_2, lam_ser_3, lam_ser_4, sqrt_series_13, sqrt_series_24] = ...
    lambda_series_num_newLam(nu, mu, N);

lam_1 = polyval(lam_ser_1, w);
lam_2 = polyval(lam_ser_2, w);
lam_3 = polyval(lam_ser_3, w);
lam_4 = polyval(lam_ser_4, w);

eq_1 = lam_eq(lam_1);
eq_2 = lam_eq(lam_2);
eq_3 = lam_eq(lam_3);
eq_4 = lam_eq(lam_4);

figure(2);
subplot(2,2,1);
contour(x, y, log10(abs(eq_1)), 'ShowText', 'on');
xlabel('Real(\omega)'); ylabel('Imag(\omega)');
title('Equation for \lambda_1(\omega)');
subplot(2,2,2);
contour(x, y, log10(abs(eq_2)), 'ShowText', 'on');
xlabel('Real(\omega)'); ylabel('Imag(\omega)');
title('Equation for \lambda_2(\omega)');
subplot(2,2,3);
contour(x, y, log10(abs(eq_3)), 'ShowText', 'on');
xlabel('Real(\omega)'); ylabel('Imag(\omega)');
title('Equation for \lambda_3(\omega)');
subplot(2,2,4);
contour(x, y, log10(abs(eq_4)), 'ShowText', 'on');
xlabel('Real(\omega)'); ylabel('Imag(\omega)');
title('Equation for \lambda_4(\omega)');