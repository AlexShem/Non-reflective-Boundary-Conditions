function u_tau = u_tau_function(h, D, C, u0)
mult = exp(-h / sqrt(D));
u0 = u0(:);
u0 = [u0(1) * mult; u0; u0(end)*mult];

a = (h^2 - 6*D)/(4*(h^2 + 3*D));
p = -(3*C)/(2*(h^4 + 3*D*h^2));
q = (6*C)/(h^4 + 3*D*h^2);
r = -(9*C)/(h^4 + 3*D*h^2);

N = length(u0)-2;
A = spdiags(ones(N, 1) .* [a 1 a], -1:1, N, N);
% A = eye(N) + a * diag(ones(N-1,1), -1) + a * diag(ones(N-1,1), 1);
B = spdiags(ones(N + 2, 1) .* [p q r q p], -2:2, N + 2, N + 2);
%     B = r*eye(N+2) + ...
%         q*(diag(ones(N+2-1,1), 1) + diag(ones(N+2-1,1), -1)) + ...
%         p*(diag(ones(N+2-2,1), 2) + diag(ones(N+2-2,1), -2));
B(1, :) = [];
B(end, :) = [];

% Boundary Conditions
B(1, :) = 0;
B(end, :) = 0;
A(1, 1) = 1;
A(1, 2) = exp(h/sqrt(D));
A(end, end) = 1;
A(end, end-1) = exp(h/sqrt(D));

% A = sparse(A);
% B = sparse(B);
u_tau = A \ (B*u0);
end
