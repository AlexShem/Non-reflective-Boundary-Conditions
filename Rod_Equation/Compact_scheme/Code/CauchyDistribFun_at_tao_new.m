function f = CauchyDistribFun_at_tao_new(a, b, x)
    f = zeros(length(x), 1);
    if a <= x(1)
        l_ind = 1;
    else
        l_ind = find(x >= a, 1);
    end
    if b >= x(end)
        r_ind = length(x);
    else
        r_ind = find(x >= b, 1);
    end
    median = (x(r_ind) + x(l_ind)) / 2;
    mlt = (x(r_ind) - x(l_ind)) / 2;
    %x = meadian + mlt * t ===> t = (x - median) / mlt
    for i = l_ind : r_ind
        cur_x = (x(i) - median) / mlt;
        f(i) = exp(-(cur_x - 1)^-2 - (cur_x + 1)^-2);
    end
    supF = max(f);
    f = f / supF;
end