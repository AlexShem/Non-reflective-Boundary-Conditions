[nu, mu] = meshgrid(linspace(1e-4, 2, 501), linspace(1e-4, 2.5, 501));

fval = zeros(size(nu));
for i = 1 : size(nu, 1)
    for j = 1 : size(mu, 1)
        a = 1+2*mu(i,j) + 3*nu(i,j); b = -mu(i,j) - 2*nu(i,j); g = 2*mu(i,j); d = -2 - 4*nu(i,j); s = nu(i,j)/2;
        fval(i,j) = b^2 - 4*s*(a-2*s);
    end
end

figure(4);
contour(nu, mu, fval, 'ShowText', 'on');
xlabel('\nu'); ylabel('\mu');
title('\beta^2 - 4\sigma (\alpha - 2\sigma)');