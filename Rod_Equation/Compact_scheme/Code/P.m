function Legandre = P(n, x)

    if n == 0
        Legandre = 1;
        return;
    end
    t = acos(x);
    extra_mul = 1;
    if ~mod(n, 2)
        c = 0;
        extra_mul = .5;
    else
        c = 1;
    end
    
    cf = ones(1, ceil((n + 1)/2));
    main_cf = 1;
    for j = 2 : n
        main_cf = main_cf * (2*j - 1) / j / 2;
    end
    
    ind = 2;
    for j = n - 2 : -2 : c
        for k = 1 : (n - j) / 2
            cf(ind) = cf(ind) * (2*k - 1) / k * (n - k + 1) / (2*n - 2*k + 1);
        end
        ind = ind + 1;
    end
    cf(end) = cf(end) * extra_mul;
    cf = cf * main_cf;
    arg = t * (n:-2:c);
    Legandre = cf * cos(arg');
end
