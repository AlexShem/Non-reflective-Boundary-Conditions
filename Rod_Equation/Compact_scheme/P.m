function Legandre = P(n, x)
%     global Pcoef;
%     if isempty(Pcoef)
%         load Pcoef.mat
%     end
%     if length(Pcoef) < n + 1 || isempty(Pcoef{n + 1})
%         for k = length(Pcoef) + 1 : n + 1
%             n1 = k - 2;
%             Pcoef{k} = [Pcoef{k - 1} 0]*(2*n1 + 1)/(n1 + 1) -  [0 0 Pcoef{k - 2}]*n1/(n1 + 1);
%         end
%         save Pcoef.mat Pcoef
%     end
%     Legandre = polyval(Pcoef{n + 1}, x);

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
%     main_cf = main_cf / 2^(n - 1);
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