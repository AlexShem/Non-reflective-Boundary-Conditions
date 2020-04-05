function I = int_sum_new(U, dx, dt, c)
    if ~ismatrix(U)
        I = NaN;
        return
    end
%     U_time = size(U, 1);
%     U_len = size(U, 2);
    f = wave_full_energy2(U, dx, dt, c);
    I = sum(f(:,2:end)+f(:,1:end-1),2)*dx*0.5;
%     I = zeros(U_time, 1);
%     for n = 1 : U_time
%         sum1 = 0;
%         for i = 1 : U_len - 1
%             f_1 = wave_full_energy(U, dx, dt, c, n, i);
%             f_2 = wave_full_energy(U, dx, dt, c, n, i + 1);
%             sum1 = sum1 + 0.5 * (f_2 + f_1) * dx;
%         end
%         I(n) = sum1;
%     end
%     f;
end