function u_tau = u_tau_function(h, D, C, f)
    mult = exp(-h / sqrt(D));
    f = f(:);
    f = [f(1) * mult; f; f(end)*mult];
    
%     a = -(- 17*h^2 + 90*D)/(56*h^2 + 180*D);
%     p = -C * (2*a + 1) / h^4;
%     q = 4*C * (2*a + 1) / h^4;
%     r = -6*C * (2*a + 1) / h^4;
    
    a = -(- h^2 + 6*D)/(4*h^2 + 12*D);
    p = -(3*C)/(2*h^2*(h^2 + 3*D));
    q = (6*C)/(h^2*(h^2 + 3*D));
    r = -(9*C)/(h^2*(h^2 + 3*D));
    
    N = length(f)-2;
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
    
%     o1 = ones(N,1);
%     A = spdiags(o1.*[a 1 a], -1:1, N, N);
%     A(1, 1) = 1; A(1, 2) = exp(h/sqrt(D));
%     A(end, end) = 1; A(end, end-1) = exp(h/sqrt(D));
% %     o1 = ones(N-2,1);
%     B = spdiags(o1.*[p q r q p], -1:3, N, N+2);
%     B(1,:) = 0; B(end, :) = 0;

    u_tau = A \ (B*f);
end