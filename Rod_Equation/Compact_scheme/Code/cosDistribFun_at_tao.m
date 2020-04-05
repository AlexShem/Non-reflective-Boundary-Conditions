function f = cosDistribFun_at_tao(a, b, x)
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
    median = (x(l_ind) + x(r_ind)) / 2;
    mlt = (x(r_ind) - x(l_ind)) / (2 * pi);
    %x = meadian + mlt * t ===> t = (x - median) / mlt, -pi < t < pi
    for i = l_ind : r_ind
        cur_x = (x(i) - median) / mlt;
        f(i) = (cos(cur_x) + 1) / 2;
    end
end