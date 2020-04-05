nu = .4; mu = .5;
% nu = 1.4; mu = .5;
% nu = .01; mu = nu + 1/12 + .05;
a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*nu; s = nu/2;

[x, y] = meshgrid(linspace(-1, 1, 81), linspace(-1, 1, 81));
w = x+1i*y;

w = linspace(1e-8, .01, 1001);

w = .01;

eq_eta = @(eta) s*(1+w.^2).*eta.^2 + (b*(1+w.^2)+g*w).*eta + a*(1+w.^2) + d*w - 2*s*(1+w.^2);

% eta_1 = -(b*(1+w.^2)+g*w)./(2*s*(1+w.^2)) - .5 * sqrt(((b*(1+w.^2)+g*w)./(s*(1+w.^2))).^2 - 4*((a*(1+w.^2)+d*w-2*s*(1+w.^2))./(s*(1+w.^2))));
% eta_2 = -(b*(1+w.^2)+g*w)./(2*s*(1+w.^2)) + .5 * sqrt(((b*(1+w.^2)+g*w)./(s*(1+w.^2))).^2 - 4*((a*(1+w.^2)+d*w-2*s*(1+w.^2))./(s*(1+w.^2))));

% eta_1 = 1./(2*s*(1+w.^2)) .* (-b*(1+w.^2) - g*w - sqrt((b^2 - 4*s*(a-2*s))*w.^4 + (2*b*g-4*d*s)*w.^3 ...
%     +(2*b^2+g^2-8*s*(a-2*s))*w.^2 + (2*b*g-4*d*s)*w + (b^2 - 4*s*(a-2*s))));
% eta_2 = 1./(2*s*(1+w.^2)) .* (-b*(1+w.^2) - g*w + sqrt((b^2 - 4*s*(a-2*s))*w.^4 + (2*b*g-4*d*s)*w.^3 ...
%     +(2*b^2+g^2-8*s*(a-2*s))*w.^2 + (2*b*g-4*d*s)*w + (b^2 - 4*s*(a-2*s))));

eps_1 = 1/(b^2-4*s*(a-2*s)) * (-b*g+2*d*s + 2*sqrt(-b*g*d*s+d^2*s^2+a*g^2*s-2*g^2*s^2));
eps_2 = 1/(b^2-4*s*(a-2*s)) * (-b*g+2*d*s - 2*sqrt(-b*g*d*s+d^2*s^2+a*g^2*s-2*g^2*s^2));

w_1 = .5*(eps_1 - sqrt(eps_1^2 - 4));
w_2 = .5*(eps_1 + sqrt(eps_1^2 - 4));
w_3 = .5*(eps_2 - sqrt(eps_2^2 - 4));
w_4 = .5*(eps_2 + sqrt(eps_2^2 - 4));

% eta_1 = 1./(2*s*(1+w.^2)) .* (-b*(1+w.^2) - g*w - sqrt((b^2-4*s*(a-2*s))*(w.^2 - eps_1*w + 1) .* (w.^2 - eps_2*w + 1)));
% eta_2 = 1./(2*s*(1+w.^2)) .* (-b*(1+w.^2) - g*w + sqrt((b^2-4*s*(a-2*s))*(w.^2 - eps_1*w + 1) .* (w.^2 - eps_2*w + 1)));

% eta_1 = 1./(2*s*(1+w.^2)) .* (-b*(1+w.^2) - g*w - sqrt(b^2-4*s*(a-2*s)) * sqrt(w.^2 - eps_1*w + 1) .* sqrt(w.^2 - eps_2*w + 1));
% eta_2 = 1./(2*s*(1+w.^2)) .* (-b*(1+w.^2) - g*w + sqrt(b^2-4*s*(a-2*s)) * sqrt(w.^2 - eps_1*w + 1) .* sqrt(w.^2 - eps_2*w + 1));

geom_sum = 0;
sum_1 = 0;
sum_2 = 0;
N = 20;
for n = 0 : N
    geom_sum = geom_sum + (-1)^n * w.^(2*n);
    sum_1 = sum_1 + P(n, eps_1/2) * w.^n;
    sum_2 = sum_2 + P(n, eps_2/2) * w.^n;
end

% eta_1 = 1./(2*s) * geom_sum .* (-b*(1+w.^2) - g*w - sqrt(b^2-4*s*(a-2*s)) * sqrt(w.^2 - eps_1*w + 1) .* sqrt(w.^2 - eps_2*w + 1));
% eta_2 = 1./(2*s) * geom_sum .* (-b*(1+w.^2) - g*w + sqrt(b^2-4*s*(a-2*s)) * sqrt(w.^2 - eps_1*w + 1) .* sqrt(w.^2 - eps_2*w + 1));

eta_1 = 1/(2*s) * geom_sum .* (-b*(1+w.^2) - g*w - sqrt(b^2 - 4*s*(a-2*s)) * (w.^2 - eps_1*w + 1) ...
    .* (w.^2 - eps_2*w + 1) .* sum_1 .* sum_2);
eta_2 = 1/(2*s) * geom_sum .* (-b*(1+w.^2) - g*w + sqrt(b^2 - 4*s*(a-2*s)) * (w.^2 - eps_1*w + 1) ...
    .* (w.^2 - eps_2*w + 1) .* sum_1 .* sum_2);

% figure(2);
% contour(x, y, log10(abs(eq_eta(eta_1))), 'ShowText', 'on');
% xlabel('Re(\omega)'); ylabel('Im(\omega)');
% title('eq for \eta_1');
%  
% figure(3);
% contour(x, y, log10(abs(eq_eta(eta_2))), 'ShowText', 'on');
% xlabel('Re(\omega)'); ylabel('Im(\omega)');
% title('eq for \eta_2');