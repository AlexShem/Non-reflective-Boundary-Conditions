[nu, mu] = meshgrid(linspace(1e-4, 2, 301), linspace(1e-4, 2.5, 301));

min_root = zeros(size(nu));
for i = 1 : size(nu, 1)
    for j = 1 : size(mu, 1)
        a = 1+2*mu(i,j) + 3*nu(i,j); b = -mu(i,j) - 2*nu(i,j); g = 2*mu(i,j); d = -2 - 4*nu(i,j); s = nu(i,j)/2;
        eq = @(w) (b^2 - 4*s*(a-2*s))*w.^4 + (2*b*g-4*d*s)*w.^3 ...
            +(2*b^2+g^2-8*s*(a-2*s))*w.^2 + (2*b*g-4*d*s)*w + (b^2 - 4*s*(a-2*s));
%         min_root(i,j) = fzero(eq, 0);
        eps_1 = 1/(b^2-4*s*(a-2*s)) * (-b*g+2*d*s+2*sqrt(-b*g*d*s+d^2*s^2+a*g^2*s-2*g^2*s^2));
        eps_2 = 1/(b^2-4*s*(a-2*s)) * (-b*g+2*d*s-2*sqrt(-b*g*d*s+d^2*s^2+a*g^2*s-2*g^2*s^2));
        w_1 = .5*(eps_1 - sqrt(eps_1^2 - 4));
        w_2 = .5*(eps_1 + sqrt(eps_1^2 - 4));
        w_3 = .5*(eps_2 - sqrt(eps_2^2 - 4));
        w_4 = .5*(eps_2 + sqrt(eps_2^2 - 4));
        w = [w_1; w_2; w_3; w_4];
        min_root(i,j) = min(abs(w));        
    end
end

figure(3);
% contour(nu, mu, min_root, 'ShowText', 'on');
contour(nu, mu, min_root, [1, 1], 'ShowText', 'on');
xlabel('\nu'); ylabel('\mu');
title('min |w_{root}|');