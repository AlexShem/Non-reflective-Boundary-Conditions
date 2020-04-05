function [y_1, y_2] = mysqrt(x)
    r = abs(x);
    phi = angle(x);
    ind = phi < 0;
    phi(ind) = phi(ind) + 2*pi;
    y_1 = sqrt(r) .* exp(1i * phi/2);
    y_2 = -y_1;
end