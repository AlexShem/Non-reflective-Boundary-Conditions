addpath('..\Code')

syms nu mu w
assume(nu > 0 & mu > 0)

a = 3/(12*mu + 4) - .5;
b = (9*nu)/(3*mu + 1) - 2;
c = 1 - (12*nu + 3)/(6*mu + 2);
d = (3*nu)/(6*mu + 2);

% eta_min = -(a*(1+w^2) + c*w) / (2*d*w) - ...
%     .5*sqrt(((a*(1+w^2) + c*w) / (d*w))^2 - 4 * (w^2 + (b-2*d)*w + 1)/(d*w));
eta_min = (-(a*(1+w^2) + c*w) - sqrt((a*(1+w^2) + c*w)^2 - 4*d*w*(w^2 + (b-2*d)*w + 1))) / (2*d*w);
eta_pl = -(a*(1+w^2) + c*w) / (2*d*w) + ...
    .5*sqrt(((a*(1+w^2) + c*w) / (d*w))^2 - 4 * (w^2 + (b-2*d)*w + 1)/(d*w));

eq_eta = @(eta) d*w*eta^2 + (a*(1+w^2) + c*w)*eta + (w^2 + (b-2*d)*w + 1);

% Chech eta_min
eq_eta_min_val = eq_eta(eta_min);
simplify(eq_eta_min_val)
% simplify(subs(eq_eta_min_val, w, 0.1))

% Chech eta_pl
eq_eta_pl_val = eq_eta(eta_pl);
simplify(eq_eta_pl_val)