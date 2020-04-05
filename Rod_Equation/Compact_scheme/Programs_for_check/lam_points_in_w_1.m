nu = .6601; mu = .2401;
a = 3/(12*mu + 4) - .5; b = 9*nu/(3*mu + 1) - 2; c = 1 - (12*nu + 3) / (6*mu + 2); d = 3*nu / (6*mu + 2);

maxR = 1;
R = linspace(1e-4, maxR, 121);
Phi = linspace(0, 2*pi, 51);
[r, phi] = meshgrid(R, Phi);

w = r .* exp(1i * phi);
[x, y] = pol2cart(phi, r);

[lam_1, lam_2, lam_3, lam_4] = lambda_series(x+1i*y, nu, mu);

plot_at_w1 = 1;

figure(2);
for i = 1 : length(Phi)
    for j = 1 : length(R)
        lam = [lam_1(i, j), lam_2(i, j), lam_3(i, j), lam_4(i, j)];
        plot(real(lam), imag(lam), 'ok'); hold on;
        
        if plot_at_w1
            plot(cos(Phi), sin(Phi), '-k');
            axis([-3 3 -3 3]);
            for k = 1 : 4
                text(real(lam(k))-.1, imag(lam(k)) - .2, ['lam_', num2str(k)]);
            end
        else
            plot(maxR*cos(Phi), maxR*sin(Phi), '-k');
            plot(0.156*cos(Phi), 0.156*sin(Phi), '-.k');
            axis([-3*maxR 3*maxR -3*maxR 3*maxR]);
            for k = 1 : 4
                text(real(lam(k))-.1*maxR, imag(lam(k)) - .2*maxR, ['lam_', num2str(k)]);
            end
        end
        plot(x(i, j), y(i, j), '*r');
        hold off;
        title({['\nu = ', num2str(nu), ', \mu = ', num2str(mu)],...
            ['\omega = ', num2str(w(i, j))]});
        pause(.1);
    end
end