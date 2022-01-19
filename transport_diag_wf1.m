clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two coupled transport equations  %
% Approx equation: M y_t + K y = 0 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  nodes:          1   2   3   4   5     2e-1 2e  2e+1               Nx-1  Nx  %
%                  |---o---|---o---|  ...  |---o---|---o---|---o---|---o---|   %
%  elements:           1       2               e              Ne-1     Ne      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% to be chosen
Ne = 100;        % number of elements
Nx = 2*Ne + 1;   % number of nodes
Nt = 2000;       % number of time instances - 1

%% parameters of the model
lam = 1;      % parameter lambda
ell = 1;      % length of the space interval
T = 2;        % end of the time interval
kappa = 1;    % coef in boundary conditions
Ni = 2;       % number of PDEs

%% we set some dependent variables
he = ell/Ne;                 % length of one element
x = linspace(0,ell,Nx);      % spatial grid with the node positions
Ntot = 2*Nx;                 % number of unknowns without BC
Nf = Ntot-1;                 % degree of freedom
ht = T/(Nt-1);               % time step
t = linspace(0, T, Nt);      % time instances

%% node numbers
NNB = reshape(1:Ntot, Ni, Nx);     % NNB(i, k) = 2*(k-1)+i

% Me = 1/30*[4, 2, -1; 2, 16, 2; -1, 2, 4];   % Element mass matrix
% Ke = 1/6*[-3, -4, 1; 4, 0, -4; -1, 4, 3];   % Element stiffness matrix

%% Element matrices
syms XI;
Ntild =  [(1-XI)*(1-2*XI), 4*XI*(1-XI), XI*(2*XI - 1)];
Dxi_Ntild = diff(Ntild, XI);

Me_syms = int((Ntild')*Ntild, 0, 1);
Ke_syms = int((Dxi_Ntild')*Ntild, 0, 1);
P1e_syms = int(Ntild(1)*(Ntild')*Ntild, 0, 1);
P2e_syms = int(Ntild(2)*(Ntild')*Ntild, 0, 1);
P3e_syms = int(Ntild(3)*(Ntild')*Ntild, 0, 1);
LL1e_syms = int((Dxi_Ntild')*Dxi_Ntild, 0, 1);

Me = double(Me_syms);
Ke = double(Ke_syms);
P1e = double(P1e_syms);
P2e = double(P2e_syms);
P3e = double(P3e_syms);
LL1e = double(LL1e_syms);

%% Assemble the matrices
M = sparse(Ntot,Ntot);    % Initialize zero matrices
K = sparse(Ntot,Ntot);

for ii = 1:Ni
    for ee = 1:Ne    
        idx = [NNB(ii, 2*ee-1), NNB(ii, 2*ee), NNB(ii, 2*ee+1)];
        M(idx, idx) = M(idx, idx) + he*Me;
        K(idx, idx) = K(idx, idx) + (-1)^ii*lam*Ke;
    end
end
% Apply the Feedback BC
K(NNB(1, Nx), NNB(1, Nx)) = K(NNB(1, Nx), NNB(1, Nx)) + lam;
K(NNB(2, Nx), NNB(1, Nx)) = K(NNB(2, Nx), NNB(1, Nx)) - kappa*lam;
K(NNB(2, 1), NNB(2, 1)) = K(NNB(2, 1), NNB(2, 1)) + lam;

% Apply the dirichlet BCs
dof = 2:Ntot;                              % degrees of fredom
M_exp = M(dof, dof); K_exp = K(dof, dof);  % explicit matrices

%% Initial data for the first equation
x0 = 0.5*ell;   % center
a = 0.2*ell;    % width
c0 = 0.1;       % magnitude
f = @(x) c0 * exp( -1/a^2*((x-x0).^2) ) ;
U0 = zeros(Ntot, 1);
for kk = 1 : Nx
    U0(NNB(1, kk), 1) = f(x(kk));
end
%plot(x, state_0(NNB(1, :), 1));

U_exp = zeros(Nf, Nt);     % Initialisation
U_exp(:, 1) = U0(dof, 1);  % set the initial data at time t1

%% Solving the ODE via implicit midpoint rule: 
M1 = M_exp + 1/2 * ht * K_exp; 
M2 = M_exp - 1/2 * ht * K_exp; 
[LL,UU,PP,QQ] = lu(M1); % LU decomposition of sparse matrix M %
                        % i.e. PMQ=LU where P and Q are unitary
                        % Thus, Mx = y equiv x = Q (U\(L\(Py)))
tic
for kk=1:Nt-1
    U_exp(:, kk+1) = QQ * ( UU\( LL\( PP * M2 * U_exp(:, kk) ) )); 
end
toc

U = zeros(Ntot, Nt);
U(2:Ntot, :) = U_exp(:, :);

tic
%% Plot

subplot(1, 2, 1);
fullState1 = zeros(Nx, Nt); fullState1(2:Nx, :) = U_exp(1:Nx-1, :);
s = surf(t,x,U(NNB(1, :), :));
s.EdgeColor = 'none';
xlabel('time'); ylabel('space'); zlabel('solution 1');
title('transport first equation');


subplot(1, 2, 2);
fullState2 = U_exp(Nx:2*Nx-1, :);
s = surf(t,x,U(NNB(2, :), :));
s.EdgeColor = 'none';
xlabel('time'); ylabel('space'); zlabel('solution 2');
title('transport second equation');
toc
