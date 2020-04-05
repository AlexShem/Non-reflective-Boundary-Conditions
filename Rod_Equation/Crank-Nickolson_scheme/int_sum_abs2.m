function I = int_sum_abs2(U, dx)
    if ~ismatrix(U)
        I = NaN;
        return
    end
    I = sqrt(sum(U(:, 2 : end - 1).^2, 2) * dx + (U(:, 1).^2 + U(:, end).^2) * 0.5 * dx);
%     I = sum(U(:, 2 : end - 1), 2) * dx + (U(:, 1) + U(:, end)) * 0.5 * dx;
end