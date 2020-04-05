H = 1;  %Initial length of the rod
L = H*10;   % Take larger domatin
h = .02;    % Step
x = -L/2 : h : L/2;

[u_0, sigma] = GaussianDistrib(-.3, .3, x);
f = GaussianDistrib(-.3, .3, [x(1)-h, x, x(end)+h]);
u_0 = u_0 .* x;
u_0 = u_0(:);
f = f .* [x(1)-h, x, x(end)+h];
f = f(:);

par = [7860, 1e-3, 210e9];  % %For Steel 
D = par(2)^2;   C = par(3) * par(2)^2 / par(1);

a = -(- 17*h^2 + 90*D)/(56*h^2 + 180*D);
p = -C * (2*a + 1) / h^4;
q = 4*C * (2*a + 1) / h^4;
r = -6*C * (2*a + 1) / h^4;

% syms a p q r;
% p = -C * (2*a + 1) / h^4;
% q = 4*C * (2*a + 1) / h^4;
% r = -6*C * (2*a + 1) / h^4;
% eq = (2*a + 1)*(-720*D*C) + 2*a*h^2 * (-360*C) - 2*p*(2*h)^6 + 2*q*h^6;
% a = solve(eq, a);

exp_as = exp(-L/2 / sqrt(D));
exp_g = (L/2)^5 * exp(-(L/2)^2 / (2*sigma));

ln_exp_g = 5*log(L/2) - (L/2)^5 / (2*sigma);
ln_exp_as = -L/2 / sqrt(D);

N = length(x);
A = eye(N) + a * diag(ones(N-1,1), -1) + a * diag(ones(N-1,1), 1);
B = r*eye(N+2) + ...
    q*(diag(ones(N+2-1,1), 1) + diag(ones(N+2-1,1), -1)) + ...
    p*(diag(ones(N+2-2,1), 2) + diag(ones(N+2-2,1), -2));
B(1,:) = []; B(end,:) = [];

% Boundary Conditions
B(1,:) = 0; B(end, :) = 0;
A(1, 1) = 1; A(1, 2) = exp(h/sqrt(D));
A(end, end) = 1; A(end, end-1) = exp(h/sqrt(D));

A = sparse(A);
B = sparse(B);
u = A \ (B*f);