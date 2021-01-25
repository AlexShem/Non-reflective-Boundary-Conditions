addpath('..\Code')

x = linspace(-.1, .1, 1001);
y = linspace(-.1, .1, 1001);
[x, y] = meshgrid(x, y);
w = x + 1i*y;
e = 600;

val = (w.^2 - 2*e*w + 1).^-.5;

N = 100;
approx = 1;
for n = 1 : N
    approx = approx + P(n, e) * w.^n;
end

df = val - approx;

figure(2)
contour(x, y, log10(abs(df)), 20, 'ShowText', 'on')
xlabel('\Re \omega');
ylabel('\Im \omega');
title(['N = ' num2str(N)]);