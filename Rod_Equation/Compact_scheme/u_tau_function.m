function u_tau = u_tau_function(h, D, C, f)
    mult = exp(-h / sqrt(D));
    f = f(:);
    f = [f(1) * mult; f; f(end)*mult];
    
    a = (h^2 - 6*D)/(4*h^2 + 12*D);
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
    u_tau = A \ (B*f);
end