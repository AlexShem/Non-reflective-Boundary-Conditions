[nu, mu] = meshgrid(linspace(1e-1, 2, 501), linspace(1e-1, 2.5, 501));
a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*nu; s = nu/2;

theta = 1./(2*s) .* (-b - sqrt(b.^2 - 4*s.*(a - 2*s)));

num = theta.^2;
denum = theta.^2 - 4;

% num = 2*theta;
% denum = sqrt(theta.^2 - 4);

d = abs(num) < abs(denum);
frac = abs(num) ./ abs(denum);

figure(3);
hold on;
contour(nu, mu, frac, 'ShowText', 'on');
xlabel('\nu'); ylabel('\mu');
title('|\theta^2| / |\theta^2 - 4|');

% contour(nu, mu, frac, [1 1], 'ShowText', 'on');