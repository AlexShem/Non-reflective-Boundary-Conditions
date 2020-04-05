function PadeErmit = Pade_Ermit_lambda_1(nu, degrees)
    [x, y] = meshgrid(linspace(-.5, .5, 100), linspace(-.5, .5, 100));
%     [x, y] = meshgrid(linspace(0.5, 150, 150), linspace(0, 150, 150));
    omega = complex(x, y);
%     omega = linspace(1e-6, 1, 10000);
%     Theta = @(omega) 1 + (1 - omega).^2 ./ (2 * nu^2 * omega);
    
    lambda_true = lam_1_true_fun(nu, omega);
%     polyCoeff = PadeCoef_Modif_const(nu, degrees);
    polyCoeff = PadeCoef_Test_const(nu, degrees);
    
    cf = zeros(max(degrees) + 1, length(degrees));
    for i = 1 : length(degrees)
        cf(1 : (degrees(i) + 1), i) = 1;
    end
    cf(cf == 1) = polyCoeff;
    
    PadeErmit = 0;
    for i = 0 : length(degrees) - 1
        PadeErmit = PadeErmit + lambda_true.^i .* polyval(flip(cf(:, length(degrees) - i)), omega);
    end
%     
%     figure(1);
%     hold on;
%     plot(log10(omega), log10(abs(PadeErmit)));
% %     plot((omega), (abs(PadeErmit)));
%     xlabel('lg(\omega)');
%     title({'lg(|Pade-Ermit|) for 4 - 3 - 1', ...
%         ['\nu = ', num2str(nu)]});

    figure(2);
%     surf(real(omega), imag(omega), abs(PadeErmit));
    contour(real(omega), imag(omega), log10(abs(PadeErmit)), 'ShowText', 'on');
    xlabel('Re(\omega)');
    ylabel('Im(\omega)');
    title({'log_{10}|Pade-Ermit| for 4 - 3 - 1', ...
        ['\nu = ', num2str(nu), '; with w = 1']});
end