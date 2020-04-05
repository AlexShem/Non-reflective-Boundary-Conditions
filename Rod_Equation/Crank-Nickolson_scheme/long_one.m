function sol = long_one()
%     clear

    n = 801;
    circ = 20;
    h = circ/(n-1);
    T = 5;

    rho = 0.5;
    R = 0.02;
%     R = 0.001;
    E = 5;

    D = R^2;
    C = E*R^2/rho;

%     nu = 0.01;
    nu = 0.256;
%     nu = 0.016;
    mu = D*h^-2;
    tau = sqrt(nu*h^4/C);
    
    x = -circ/2:h:circ/2;

    sol = zeros(floor(T/tau), length(x));
    
    a = 3/(12*mu + 4) - 1/2;
    b = (9*nu)/(3*mu + 1) - 2;
    c = 1 - (12*nu + 3)/(6*mu + 2);
    d = (3*nu)/(6*mu + 2);

    p0 = 1/48 - (3*nu - 3)/(48*(3*mu + 1));
    q0 = (2*nu + 5/3)/(48*(3*mu + 1)) - 1/72;
    r0 = 1/288 - (nu + 1/3)/(96*(3*mu + 1));
    p1 = 5/24 - (15*nu - 15)/(24*(3*mu + 1));
    q1 = (10*nu + 25/3)/(24*(3*mu + 1)) - 5/36;
    r1 = 5/144 - (5*nu + 5/3)/(48*(3*mu + 1)); 

    z = 0;
    a02 = -(15*mu + 25*nu + 36*mu*z + 60*nu*z - 108*mu^2*z + 180*mu*nu*z)/(10*(6*mu + 5*nu));
    a03 = (10*mu + 25*nu + 24*mu*z + 60*nu*z + 360*mu*nu*z)/(45*(6*mu + 5*nu));
    a11 = -(4*(15*mu + 5*nu - 18*nu*z - 162*mu*nu*z))/(5*(6*mu + 5*nu));
    a12 = (15*mu + 10*nu + 36*mu*z + 24*nu*z - 108*mu^2*z - 144*mu*nu*z)/(5*(6*mu + 5*nu));
    a13 = -(4*(5*mu + 5*nu + 12*mu*z + 12*nu*z + 18*mu*nu*z))/(45*(6*mu + 5*nu));
    p00 = z/20 + (z/10 - (3*nu*z)/4 + 1/24)/(6*mu + 5*nu);
    p01 = z/10;
    p10 = (z + 9*mu*z + 5/12)/(6*mu + 5*nu) - z;
    p11 = z;

    % 0 for values and derivatives on segment ends
    U1 = diag(a*ones(n-1,1),-1) + diag(a*ones(n-1,1),1) + diag(ones(n,1));
    U1(1, :) = 0; U1(2, :) = 0; U1(1,1) = 1; U1(2, 2) = 1; U1(2,3) = a02; U1(2,4) = a03;
    U1(end, :) = 0; U1(end-1, :) = 0; U1(end,end) = 1; U1(end-1, end-1) = 1; U1(end-1,end-2) = a02; U1(end-1,end-3) = a03;
    U2 = diag(d*ones(n-2,1),-2) + diag(d*ones(n-2,1),2) + ...
        diag(c*ones(n-1,1),-1) + diag(c*ones(n-1,1),1) + ...
        diag(b*ones(n,1));

    U2(1, :) = 0; U2(2, :) = 0; U2(2,2) = a11; U2(2,3) = a12; U2(2,4) = a13;  
    U2(end, :) = 0; U2(end-1, :) = 0; U2(end-1, end-1) = a11; U2(end-1,end-2) = a12; U2(end-1,end-3) = a13;

    F1 = diag(r0*ones(n-2,1),-2) + diag(r0*ones(n-2,1),2) + ...
        diag(q0*ones(n-1,1),-1) + diag(q0*ones(n-1,1),1) + ...
        diag(p0*ones(n,1));
    F1(1,:) = 0; F1(2,:) = 0; 
    F1(2,1) = p00; F1(2,2) = p01;
    F1(end,:) = 0; F1(end-1,:) = 0; 
    F1(end-1,end) = p00; F1(end-1,end-1) = p01;

    F2 = diag(r1*ones(n-2,1),-2) + diag(r1*ones(n-2,1),2) + ...
        diag(q1*ones(n-1,1),-1) + diag(q1*ones(n-1,1),1) + ...
        diag(p1*ones(n,1));
    F2(1,:) = 0; F2(2,:) = 0; 
    F2(2,1) = p10; F2(2,2) = p11;
    F2(end,:) = 0; F2(end-1,:) = 0; 
    F2(end-1,end) = p10; F2(end-1,end-1) = p11;

    F1 = F1*tau^2; F2 = F2*tau^2;

    % % periodical BC
    % U1 = diag(a*ones(n-1,1),-1) + diag(a*ones(n-1,1),1) + diag(ones(n,1));
    % U1(end,1) = a; U1(1,end) = a;
    % 
    % U2 = diag(d*ones(n-2,1),-2) + diag(d*ones(n-2,1),2) + ...
    %     diag(c*ones(n-1,1),-1) + diag(c*ones(n-1,1),1) + ...
    %     diag(b*ones(n,1));
    % U2(end,1) = c; U2(1,end) = c;
    % U2(end,2) = d; U2(2,end) = d;
    % U2(end-1,1) = d; U2(1,end-1) = d;
    % 
    % F1 = diag(r0*ones(n-2,1),-2) + diag(r0*ones(n-2,1),2) + ...
    %     diag(q0*ones(n-1,1),-1) + diag(q0*ones(n-1,1),1) + ...
    %     diag(p0*ones(n,1));
    % F1(end,1) = q0; F1(1,end) = q0;
    % F1(end,2) = r0; F1(2,end) = r0;
    % F1(end-1,1) = r0; F1(1,end-1) = r0;
    % F1 = F1*tau^2;
    % 
    % F2 = diag(r1*ones(n-2,1),-2) + diag(r1*ones(n-2,1),2) + ...
    %     diag(q1*ones(n-1,1),-1) + diag(q1*ones(n-1,1),1) + ...
    %     diag(p1*ones(n,1));
    % F2(end,1) = q1; F2(1,end) = q1;
    % F2(end,2) = r1; F2(2,end) = r1;
    % F2(end-1,1) = r1; F2(1,end-1) = r1;
    % F2 = F2*tau^2;


    u2 = x*0;
    for cnt = 1:length(x)
        if abs(x(cnt)) <= pi/2
            u2(cnt) = cos(x(cnt))^2;
        end
    end

    sol(1:2, :) = repmat(u2, 2, 1);
    
    t = tau;
    u3 = u2.';
    uu = u2.';
    u1 = u2.';
    u2 = u2.';
    % U1 = sparse (U1);
    U1i = U1^-1;

%     gfx = plot(x,u3,'LineWidth',2);
    % tic
    for cnt = 3:T/tau    
%         set(gfx, 'YData', u3);
%     %     plot(x,u3,'LineWidth',2)
%         title(t)
%         grid on
%         drawnow
%         axis([-circ/2 circ/2 -0.3 1.2])
%     %     pause(1e-4); 

        f3 = x.'*0;
        f2 = x.'*0;
        f1 = x.'*0;

        right_hand = F1*(f1+f3) + F2*f2 - U2*u2;
    %     u3 = U1 \ right_hand - u1;
        u3 = U1i * right_hand - u1;
        sol(cnt, :) = u3;

        uu = u1;
        u1 = u2;
        u2 = u3;
        t = t + tau;
    end
    % toc
end