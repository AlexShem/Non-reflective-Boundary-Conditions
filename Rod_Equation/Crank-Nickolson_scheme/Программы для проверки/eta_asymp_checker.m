[nu, mu] = meshgrid(linspace(1e-4, 2, 501), linspace(1e-4, 2.5, 501));
a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*nu; s = nu/2;

th_1 = 1./(2*s) .* (-b - sqrt(b.^2 - 4*s.*(a-2*s)));
th_2 = 1./(2*s) .* (-b + sqrt(b.^2 - 4*s.*(a-2*s)));

val_1 = .25 * th_1.^2 - 1;
val_2 = .25 * th_2.^2 - 1;

figure(1);
contour(nu, mu, abs(val_1), 'ShowText', 'on');

figure(2);
contour(nu, mu, abs(val_2), 'ShowText', 'on');

eps_1 = th_1 ./ sqrt(th_1.^2 - 4);
eps_2 = th_2 ./ sqrt(th_2.^2 - 4);