% fraction series
nu = 9000; mu = 10000;
a = 3/(12*mu + 4) - .5; b = 9*nu/(3*mu + 1) - 2; c = 1 - (12*nu + 3) / (6*mu + 2); d = 3*nu / (6*mu + 2);

[x, y] = meshgrid(linspace(-1, 1, 101), linspace(-1, 1, 101));
w = complex(x, y);

fr = 1./(w.^2 + (b-2*d)*w + 1);

A = -1/sqrt((b-2*d)^2 -4);
B = -A;
w_1 = .5 * (-b+2*d - sqrt((b-2*d)^2 - 4));
w_2 = .5 * (-b+2*d + sqrt((b-2*d)^2 - 4));

fr_ser = 0;
N = 100;
for n = 0 : N
    fr_ser = fr_ser + (-1/w_1^(n+1) + 1/w_2^(n+1)) * w.^n;
end
fr_ser = A * fr_ser;

figure(1);
contour(real(w), imag(w), log10(abs(fr - fr_ser)), 'ShowText', 'on');