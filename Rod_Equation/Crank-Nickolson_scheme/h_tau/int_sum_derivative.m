function I = int_sum_derivative(U, dx, dt, rho, R, E)
%     dtU = (U(3:end, 2:end-1) - U(1:end-2, 2:end-1))/(2*dt);
%     dxdtU = (U(3:end, 3:end) + U(1:end-2, 1:end-2) - U(3:end, 1:end-2) - U(1:end-2, 3:end))/(4*dx*dt);
%     d2xU = (U(2:end-1, 3:end) - 2*U(2:end-1, 2:end-1) + U(2:end-1, 1:end-2))/dx^2;
%     
%     K_fun = (dtU.^2 + R^2 * dxdtU.^2);
%     K = .5 * rho * (sum(K_fun(:, 2:end-1), 2) + .5*K_fun(:, 1) + .5*K_fun(:, end))*dx;
%     
%     P_fun = d2xU.^2;
%     P = .5 * E*R^2 * (sum(P_fun(:, 2:end-1), 2) + .5*P_fun(:, 1) + .5*P_fun(:, end))*dx;
%     I = K + P;

    dtU = (U(2:end, :) - U(1:end-1, :)) / dt;
%     dtU(:, 1) = []; dtU(:, end) = [];
    dxdtU = (U(2:end, 3:end) - U(1:end-1, 3:end) - U(2:end, 1:end-2) + U(1:end-1, 1:end-2))/(2*dx*dt);
    d2xU = (U(2:end, 3:end) - 2*U(2:end, 2:end-1) + U(2:end, 1:end-2) + ...
        U(1:end-1, 3:end) - 2*U(1:end-1, 2:end-1) + U(1:end-1, 1:end-2))/(2*dx^2);
    
    dxdtU_0 = (U(2:end, 2) - U(1:end-1, 2) - U(2:end, 1) + U(1:end-1, 1))/(dx*dt);
    dxdtU_N = (-U(2:end, 2) + U(1:end-1, 2) + U(2:end, 1) - U(1:end-1, 1))/(dx*dt);
    
    theta = rho * dtU(:, 2:end-1).^2 + rho*R^2 * dxdtU.^2 + E*R^2 * d2xU.^2;
    theta_0 = rho * dtU(:, 1).^2 + rho*R^2 * dxdtU_0.^2 + E*R^2 * d2xU(:, 1).^2;
    theta_N = rho * dtU(:, end).^2 + rho*R^2 * dxdtU_N.^2 + E*R^2 * d2xU(:, end).^2;
    
    I = sqrt(dx * (.5*(theta_0 + theta_N) + sum(theta, 2)));
end