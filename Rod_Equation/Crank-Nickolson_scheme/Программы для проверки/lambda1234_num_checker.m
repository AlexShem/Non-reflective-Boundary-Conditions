nu = .05; mu = .15;
N_lam = 50;
a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*nu; s = nu/2;

[x, y] = meshgrid(linspace(-1, 1, 101), linspace(-1, 1, 101));
w = x+1i*y;

[lam_ser_1, lam_ser_2, lam_ser_3, lam_ser_4] = lambda_series_num(nu, mu, N_lam);
lam = {lam_ser_1, lam_ser_2, lam_ser_3, lam_ser_4};
eq = cell(1, 4);

figure(1);
for i = 1 : 4
    lam_w = polyval(lam{i}, w);
    eq{i} = s*(1 + w.^2) .* lam_w.^4 + s*(1 + w.^2) + (b*(1+w.^2)+g*w).*lam_w.^3 + ...
        (b*(1+w.^2)+g*w).*lam_w + (a*(1+w.^2)+d*w).*lam_w.^2;
    subplot(2, 2, i);
    contour(x, y, log10(abs(eq{i})), 'ShowText', 'on');
    title(['eq for \lambda_' num2str(i)]);
    xlabel('Re(\omega)'); ylabel('Im(\omega)');
end