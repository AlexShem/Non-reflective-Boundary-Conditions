function f = wave_full_energy2(U, dx, dt, c)
    f_dx = zeros(size(U));
    f_dx(:,2:end-1) = (U(:, 3:end) - U(:, 1:end-2)) / (2 * dx);
    f_dx(:,1) = (-3 * U(:, 1) + 4 * U(:, 2) - U(:, 3)) / (2 * dx);
    f_dx(:,end)= (3 * U(:, end - 2) - 4 * U(:, end - 1) + U(:, end)) / (2 * dx);

    f_dt = zeros(size(U));
    f_dt(2:end-1,:) = (U(3:end,:) - U(1:end-2,:)) / (2 * dt);
    f_dt(1,:) = (-3 * U(1,:) + 4 * U(2,:) - U(3,:)) / (2 * dt);
    f_dt(end,:)= (3 * U(end - 2,:) - 4 * U(end - 1,:) + U(end,:)) / (2 * dt);

    f = c^2 * f_dx.^2 + f_dt.^2;
end