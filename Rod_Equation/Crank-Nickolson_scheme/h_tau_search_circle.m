warning('off','all');

%%---------- Default ----------
% [h, tau] = meshgrid(...
%     .0005 : .00025 : .03,...
%     .00001 : .00000125 : .0002);

%%---------- Precise calculations ----------
[h, tau] = meshgrid(...
    linspace(.0001, .01, 250),...
    linspace(.00001, .1, 250));
% [h, tau] = meshgrid(...
%     linspace(.0001, .07, 850),...
%     linspace(.00001, .0006, 750));

%%----------Testing the programm ----------
% [h, tau] = meshgrid(...
%     linspace(.001, 1.0, 30),...
%     linspace(.0001, 1.0, 30));
% [h, tau] = meshgrid(...
%     linspace(.0001, 1, 40),...
%     linspace(.00001, 1, 40));

par = [7860, 1e-3, 210e9];  % %For Steel
D = par(2)^2;   C = par(3) * par(2)^2 / par(1);

nx = size(h, 2); nt = size(tau, 1);
stable = zeros(nt, nx);

tic
for j = 1 : nx
    for i = 1 : nt
        nu = C * tau(i, j).^2 ./ h(i, j).^4; % Параметр Куранта
        mu = D ./ h(i, j).^2;
        a = 1 + 2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*mu; s = nu/2;
        theta = @(z) (2*g*z + d)./(2*s*(2*z.^2 - 1) + 2*b*z + a);
        
        options = optimset('TolX', 1e-8);
        %         options = optimset('PlotFcns',@optimplotfval);
        [zmin, th_min, exitflag, output] = fminbnd(theta, -1, 1, options);
        [zmax, th_max, exitflag, output] = fminbnd(@(z) -1*theta(z), -1, 1, options);
        th_max = -th_max;
        if th_min < -1 || th_max > 1
            stable(i, j) = false;
        else
            stable(i, j) = true;
        end
    end
    if j == floor(nx/2)
        disp('50% for circle');
    elseif j == floor(nx/4)
        disp('25% for circle');
    elseif j == floor(.75*nx)
        disp('75% for circle');
    elseif j == nx
        disp('100% for circle');
    end
end
toc

figure(19);
imagesc(h(1,:), tau(:,1), stable); axis xy; colormap([0 0 0; 1 1 1]);
xlabel('$h$ (m)', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('$\tau$ (s)', 'Interpreter', 'latex', 'FontSize', 15);
title('Circle CN');

filename = 'CN_circle.mat';
save(filename, 'h', 'tau', 'stable');