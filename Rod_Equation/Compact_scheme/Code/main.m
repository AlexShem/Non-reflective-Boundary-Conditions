%% Define parameters

% Length of the rod
L = 1;
%Integration time
T = .5;

% Parameters of the steel rod
rho = 7860;
R = 1e-2;
E = 210e9;
rod_par = [rho, R, E];
D = R^2;
C = E*R^2/rho;

h = 0.01;
tau = 0.0001;
% Calculate the number of spatial steps
Nx = ceil(L / h) + 1;
% Adjust step h
h = L / (Nx - 1);
% Calculate the number of temporal steps
Nt = ceil(T / tau) + 1;
% Adjust step h
tau = T / (Nt - 1);

% Calculate dimensionless parameters
nu = C * tau^2 / h^4;
mu = D / h^2;

% Calculate the coefficients of the compact scheme
a = 3/(12*mu + 4) - .5;
b = (9*nu)/(3*mu + 1) - 2;
c = 1 - (12*nu + 3)/(6*mu + 2);
d = (3*nu)/(6*mu + 2);



%% Initial Conditions

% Spatial grid
x = linspace(-L/2, L/2, Nx);
% Temporal grid
t = linspace(0, T, Nt);

% Initial condition
U_0 = GaussianDistrib(-.5, .5, x);
U_0 = x.*U_0;

U_1 = 0 * U_0;
U_2 = u_tau_function(h, D, C, U_0);
U_tau = U_0 + U_1*tau + U_2*tau^2/2;