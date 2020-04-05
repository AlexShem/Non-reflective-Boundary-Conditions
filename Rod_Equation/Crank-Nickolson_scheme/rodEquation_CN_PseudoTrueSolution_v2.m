function U = rodEquation_CN_PseudoTrueSolution_v2(steps, par, kurant, T, u_0, u_tau, domain)
%     Set the parameters
%   par = [rho, R, E]

    h = steps(1); tau = steps(2);
    D = par(2)^2; C = par(3) * par(2)^2 / par(1);

%     nu = C * tau^2 / h^4; % Параметр Куранта
%     mu = D / h^2;

    nu = kurant(1); mu = kurant(2);
    a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*mu; s = nu/2;
    
%     L = 10;  %Длина струны, м.
    
    range = domain(end)-domain(1);
%     exdend = range * 80;
    exdend = range * 140;
%     exdend = range * 260;
    u_0 = [zeros(size(domain(1):-h:(domain(1)-exdend))),...
        u_0.', ...
        zeros(size((domain(end):h:(domain(end)+exdend)))) ];
    u_tau = [zeros(size(domain(1):-h:(domain(1)-exdend))),...
        u_tau.', ...
        zeros(size((domain(end):h:(domain(end)+exdend)))) ];
        
    Nx = length(u_0); % Кол-во точек по пространству
    Nt = length(0 : tau : T);   % Кол-во точек по времени
         
    U = zeros(Nt, Nx);
    U(1, :) = u_0';
    U(2, :) = u_tau';

%     U_now = diag(a*ones(1, Nx)) + diag(b * ones(1, Nx-1), 1) + diag(b * ones(1, Nx-1), -1) + diag(s * ones(1, Nx-2), 2) + diag(s * ones(1, Nx-2), -2);
    
    U_now = spdiags(ones(Nx, 1) .* [s b a b s], -2:2, Nx, Nx);
    U_at_pTime = spdiags(ones(Nx, 1) .* [g d g], -1:1, Nx, Nx);
    U_at_ppTime = U_now;
% %  Diriclet Boundary
    U_now([1,2,end-1, end],:) = []; U_now(:,[1,2,end-1,end]) = [];
%     U_at_pTime = diag(d * ones(1, Nx)) + diag(g * ones(1, Nx-1), 1) + diag(g * ones(1, Nx-1), -1);
    U_at_pTime([1,2,end-1, end],:) = []; U_at_pTime(:,[1,2,end-1,end]) = [];
%     U_at_ppTime = diag(a*ones(1, Nx)) + diag(b * ones(1, Nx-1), 1) + diag(b * ones(1, Nx-1), -1) + diag(s * ones(1, Nx-2), 2) + diag(s * ones(1, Nx-2), -2);
    U_at_ppTime([1,2,end-1, end],:) = []; U_at_ppTime(:,[1,2,end-1,end]) = [];

%     U_now = sparse(U_now);
%     U_at_pTime = sparse(U_at_pTime);
%     U_at_ppTime = sparse(U_at_ppTime);
    for n = 3 : Nt
        rightPart = -U_at_pTime * U(n - 1, 3:end-2)' - U_at_ppTime * U(n - 2, 3:end-2)';
        U(n, 3:end-2) = U_now \ rightPart;
    end

% % Animate in programm
%     figure(6)
%     pause(1);
%     for i = 1 : 1 : Nt
%         plot(linspace(-exdend/2, exdend/2, Nx), U(i,:));
%         axis([-exdend/2, exdend/2, min(u_0), max(u_0)]);
%         title(['t = ' num2str((i-1)*tau)]);
%         pause(.01);
%     end

end