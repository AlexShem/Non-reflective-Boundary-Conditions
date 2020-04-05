warning('off','all');

par = [7860, 1e-3, 210e9];
rho = par(1); R = par(2); E = par(3);

% % No constant condition
degrees = {%[4, 4, 8, 8];
     [8, 8, 6, 6];
    [5, 3, 8, 8];
%     [5, 3, 9, 7]
    };

%%---------- Precise calculations ----------
% [h, tau] = meshgrid(...
%     .0005 : .000125 : .1,...
%     .00001 : .0000004 : .00025);
[h, tau] = meshgrid(...
    linspace(.0001, .07, 750),...
    linspace(.00001, .0006, 750));

%%----------Testing the programm ----------
% [h, tau] = meshgrid(...
%     linspace(.0001, .07, 30),...
%     linspace(.00001, .0006, 30));
% [h, tau] = meshgrid(...
%     linspace(.0001, .25, 20),...
%     linspace(.00001, .0025, 20));
% [h, tau] = meshgrid(...
%     linspace(.0001, .5, 20),...
%     linspace(.00001, .005, 20));
% h = .02; tau = 1.6e-4;

nx = size(h, 2); nt = size(tau, 1);
gd_max_H = zeros(size(h));
gd_max_L2 = zeros(size(h));
gd_max_Ch = zeros(size(h));
bc_exist = false(size(h));

for k = 1 : length(degrees)
    polynomDEG_bor = degrees{k};
    tic
%     for j = 1 : nx        
%         parfor i = 1 : nt
%             [bc_exist(i,j),gd_max_H(i,j),gd_max_L2(i,j),gd_max_Ch(i,j)]=calconetauh (tau(i,j),h(i,j),polynomDEG_bor, rho, R, E);  
%         end        
%         disp([num2str(100*j/nx) '% for ', num2str(polynomDEG_bor)]);
%         toc
%     end
    parfor i = 1 : numel(tau)
        [bc_exist(i),gd_max_H(i),gd_max_L2(i),gd_max_Ch(i)]=calconetauh (tau(i),h(i),polynomDEG_bor, rho, R, E);  
        disp([num2str(100*i/numel(tau)) '% for ', num2str(polynomDEG_bor)]);
    end        
    
    figure(3*k - 2);
    imagesc(h(1,:), tau(:,1), gd_max_H); axis xy; colormap([0 0 0; 1 1 1]);
    xlabel('$h$ (m)', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('$\tau$ (s)', 'Interpreter', 'latex', 'FontSize', 15);
    title(['H ', num2str(polynomDEG_bor)]);
    
    figure(3*k - 1);
    imagesc(h(1,:), tau(:,1), gd_max_Ch); axis xy; colormap([0 0 0; 1 1 1]);
    xlabel('$h$ (m)', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('$\tau$ (s)', 'Interpreter', 'latex', 'FontSize', 15);
    title(['C max ', num2str(polynomDEG_bor)]);
    
    figure(3*k);
    imagesc(h(1,:), tau(:,1), gd_max_L2); axis xy; colormap([0 0 0; 1 1 1]);
    xlabel('$h$ (m)', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel('$\tau$ (s)', 'Interpreter', 'latex', 'FontSize', 15);
    title(['L^2 ', num2str(polynomDEG_bor)]);
    
    filename = [num2str(polynomDEG_bor), 'no_const.mat'];
    save(filename, 'h', 'tau', 'gd_max_L2', 'gd_max_Ch', 'bc_exist');
    filename = [num2str(polynomDEG_bor), 'no_const_energy.mat'];
    save(filename, 'h', 'tau', 'gd_max_H');
end