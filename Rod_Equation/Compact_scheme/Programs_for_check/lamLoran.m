function [LegandrSum, GeomSum_1] = lamLoran(w, nu, mu, deg)
%     w_1_inf = 1 - 12*nu / (1 + 6*nu);
    
    LegandrSum = 0;
    leg_coef = 1; leg_mul = sqrt((1 - 3*nu) / (1 + 3*nu));
    GeomSum_1 = 0;
    geom_1_coef = 1; geom_1_mul = (1 + 6*nu) / (1 - 6*nu);
    for n = 0:deg
        LegandrSum = LegandrSum + P(n, mu) * leg_coef * w^n;
        GeomSum_1 = GeomSum_1 + geom_1_coef * w^n;
        leg_coef = leg_coef * leg_mul;
        geom_1_coef = geom_1_coef * geom_1_mul;
    end
end