nu = 1; mu = 1.2;
a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*nu; s = nu/2;
eq = @(eps) (b^2 - 4*s*(a-2*s)) * eps.^2 + (2*b*g - 4*d*s) * eps + g^2;

eps = linspace(-7, 7, 501);
figure(1);
plot(eps, eq(eps));

eps_1 = 1/(b^2-4*s*(a-2*s)) * (-b*g+2*d*s+2*sqrt(-b*g*d*s+d^2*s^2+a*g^2*s-2*g^2*s^2));
eps_2 = 1/(b^2-4*s*(a-2*s)) * (-b*g+2*d*s-2*sqrt(-b*g*d*s+d^2*s^2+a*g^2*s-2*g^2*s^2));

eq(eps_1)
eq(eps_1)