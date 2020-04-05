function U = rodEquation_CN_PseudoTrueSolution(steps, kurant, T, u_0, u_tau, domain)

    h = steps(1); tau = steps(2);

    nu = kurant(1); mu = kurant(2);
    a = 1+2*mu + 3*nu; b = -mu - 2*nu; g = 2*mu; d = -2 - 4*mu; s = nu/2;
    
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
        
    Nx = length(u_0);
    Nt = length(0 : tau : T);
         
    U = zeros(Nt, Nx);
    U(1, :) = u_0';
    U(2, :) = u_tau';
    
    U_now = spdiags(ones(Nx, 1) .* [s b a b s], -2:2, Nx, Nx);
    U_at_pTime = spdiags(ones(Nx, 1) .* [g d g], -1:1, Nx, Nx);
    U_at_ppTime = U_now;
% %  Diriclet Boundary
    U_now([1,2,end-1, end],:) = []; U_now(:,[1,2,end-1,end]) = [];
    U_at_pTime([1,2,end-1, end],:) = []; U_at_pTime(:,[1,2,end-1,end]) = [];
    U_at_ppTime([1,2,end-1, end],:) = []; U_at_ppTime(:,[1,2,end-1,end]) = [];

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
