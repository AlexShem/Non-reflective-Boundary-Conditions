function [eta_plus, eta_minus] = eta_true(nu, mu, w)
    a = 3/(12*mu + 4) - .5;
    b = 9*nu/(3*mu + 1) - 2;
    c = 1 - (12*nu + 3) / (6*mu + 2);
    d = 3*nu / (6*mu + 2);

    theta = ((a*(1+w.^2)+c*w) ./ (d*w)).^2 - 4 * (w.^2 + (b-2*d)*w+1) ./ (d*w);
    [rt_pos, rt_neg] = mysqrt(theta);
    pos_ind = real(rt_pos) >= 0;
    eta = -(a*(1+w.^2) + c*w) ./ (2*d*w);
    
    eta_plus = (eta + .5 * rt_pos) .* pos_ind + (eta + .5 * rt_neg) .* (1 - pos_ind);
    eta_minus = (eta - .5 * rt_pos) .* pos_ind + (eta - .5 * rt_neg) .* (1 - pos_ind);
    
%     eta_plus = eta + .5 * rt_pos;
%     eta_minus = eta - .5 * rt_pos;
end