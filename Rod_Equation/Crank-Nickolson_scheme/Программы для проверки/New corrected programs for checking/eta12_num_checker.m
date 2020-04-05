nu = .05; mu = .15;
N_eta = 80;
a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*nu; s = nu/2;

[x, y] = meshgrid(linspace(-1, 1, 101), linspace(-1, 1, 101));
w = x+1i*y;

[eta_ser_minus, eta_ser_plus] = eta_series_num(nu, mu, N_eta);
eta = {eta_ser_minus, eta_ser_plus};
eq = cell(size(eta));

figure(1);
for i = 1 : 2
    eta_w = polyval(eta{i}, w);
    eq{i} = s*(1 + w.^2) .* eta_w.^2 + (b*(1+w.^2)+g*w).*eta_w + ...
        a*(1 + w.^2) + d*w - 2*s*(1 + w.^2);
    subplot(1, 2, i);
    contour(x, y, log10(abs(eq{i})), 'ShowText', 'on');
    title(['eq for \eta_' num2str(i)]);
    xlabel('Re(\omega)'); ylabel('Im(\omega)');
end