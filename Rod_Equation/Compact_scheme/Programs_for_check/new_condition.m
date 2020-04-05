% nu = linspace(1, 5e4, 101);
% mu = linspace(1, 1e5, 101);
nu = linspace(1e-4, 1.5, 201);
mu = linspace(1e-4, 2, 201);
phi = linspace(0, 2*pi, 41);

% [nu, mu] = meshgrid(nu, mu);
a = @(n, m) 3./(12*m + 4) - .5;
b = @(n, m) 9*n./(3*m + 1) - 2;
c = @(n, m) 1 - (12*n + 3) ./ (6*m + 2);
d = @(n, m) 3*n ./ (6*m + 2);

% f = @(n, m, p) abs((2*a(n, m)*cos(p)+1) ./ (4*d(n, m)*cos(p).^2 + 2*c(n, m)*cos(p) + b(n, m)-2*d(n, m)));
f = @(n, m, p) abs((4*d(n, m)*cos(p).^2 + 2*c(n, m)*cos(p) + b(n, m)-2*d(n, m)) ./ (2*a(n, m)*cos(p)+1));
f_val_max = zeros(length(nu), length(mu));

for i = 1 : length(nu)
    for j = 1 : length(mu)
        for k = 1 : length(phi)
            f_val = f(nu(i), mu(j), phi(k));
            if f_val_max(i,j) < f_val
                f_val_max(i,j) = f_val;
            end
        end
    end
end
g = f_val_max - 2;
[nu, mu] = meshgrid(nu, mu);
figure(1);
% contour(nu, mu, g, [0,0], 'ShowText', 'on');
contour(nu, mu, g, 'ShowText', 'on');
xlabel('\nu'); ylabel('\mu');