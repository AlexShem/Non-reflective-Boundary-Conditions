function wave = wave_true_sol(x, t, c)
    % For u_0(x) = 0.5 * (cos(2*pi*x - pi) + 1);
    % and x \in [0, 1];
    
    % phi - going to the left
    % psi - going to the right
    phi = zeros(size(x));
    psi = zeros(size(x));
    x1 = x - c*t; x2 = x + c*t;
    ind_x1 = find(x1 >= 0 & x1 <= 1);
    ind_x2 = find(x2 >= 0 & x2 <= 1);
    
    phi(ind_x1) = 1/4 * (cos(2*pi*x1(ind_x1) - pi) + 1);
    psi(ind_x2) = 1/4 * (cos(2*pi*x2(ind_x2) - pi) + 1);

    wave = phi + psi;
end