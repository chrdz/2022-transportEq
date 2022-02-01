clear all
close all
clc

% "physical system" is not diagonal
% we diagonalize before discretizing and solving
% then we go back to the "physical unknown"

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
Nt = 2000*5;       % number of time instances - 1

%% parameters of the model
%lam = 1;      % parameter lambda
ell = 1;      % length of the space interval
T = 7;        % end of the time interval
Ni = 2;       % number of PDEs
massmat = eye(2);

% L = 1/sqrt(2)*[[-1, 1]; [1, 1]];
Dm = -1;
Dp = 1;
D = [[Dm, 0]; [0, Dp]];
% C1 = 1/sqrt(2)*(L')*[[-Dm, 0]; [0, Dp]]*ones(2);
% C2 = 1/sqrt(2)*(L')*D*[[-1, 1]; [-1, 1]];
% C3 = [-1; 1]*Dm;

P1 = [[-Dm, 0]; [0, Dp]]*[[0, 1]; [0, 1]];
P2 = -D*[[1, 0]; [1, 0]];
P3 = sqrt(2)*[1; 0]*Dm;

%% we set some dependent variables
he = ell/Ne;                 % length of one element
x = linspace(0,ell,Nx);      % spatial grid with the node positions
Ntot = 2*Nx;                 % number of unknowns without BC
Nf = Ntot;                 % degree of freedom
ht = T/(Nt-1);               % time step
t = linspace(0, T, Nt);      % time instances

%% node numbers
NNB = reshape(1:Ntot, Ni, Nx);     % NNB(i, k) = 2*(k-1)+i

Me = 1/30*[4, 2, -1; 2, 16, 2; -1, 2, 4];   % Element mass matrix
Ke = 1/6*[-3, -4, 1; 4, 0, -4; -1, 4, 3];   % Element stiffness matrix

%% Assemble the matrices
M = sparse(Ntot,Ntot);    % Initialize zero matrices
K = sparse(Ntot,Ntot);
W = sparse(Ntot, 1);

for ii = 1:Ni
    for jj = 1:Ni
        for ee = 1:Ne    
            idxR = [NNB(ii, 2*ee-1), NNB(ii, 2*ee), NNB(ii, 2*ee+1)];
            idxC = [NNB(jj, 2*ee-1), NNB(jj, 2*ee), NNB(jj, 2*ee+1)];
            M(idxR, idxC) = M(idxR, idxC) + he*massmat(ii, jj)*Me;
            K(idxR, idxC) = K(idxR, idxC) - D(ii, jj)*Ke;
        end
        K(NNB(ii, Nx), NNB(jj, Nx)) = K(NNB(ii, Nx), NNB(jj, Nx)) + P1(ii, jj);
        K(NNB(ii, 1), NNB(jj, 1)) = K(NNB(ii, 1), NNB(jj, 1)) + P2(ii, jj);
    end
    W(NNB(ii, Nx), 1) = W(NNB(ii, Nx), 1) + P3(ii, 1);
end
% Apply BC
% K(NNB(2, Nx), NNB(1, Nx)) = K(NNB(2, Nx), NNB(1, Nx)) + 1;
% W(NNB(1, Nx), 1) = W(NNB(1, Nx), 1) + 1;

% Apply the dirichlet BCs
dof = 1:Ntot;                              % degrees of fredom
M_exp = M(dof, dof); K_exp = K(dof, dof);  % explicit matrices
W = W(dof, :); 

%% Initial data for the first equation
U0 = zeros(Ntot, 1);
% x0 = 0.5*ell;   % center
% a = 0.2*ell;    % width
% c0 = 0.1;       % magnitude
% f = @(x) c0 * exp( -1/a^2*((x-x0).^2) ) ;
% 
% for kk = 1 : Nx
%     U0(NNB(1, kk), 1) = f(x(kk));
% end
% plot(x, U0(NNB(1, :), 1));

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
    F_k = f_ext(t(kk)+ht/2);
    U_exp(:, kk+1) = QQ * ( UU\( LL\( PP * ( M2*U_exp(:, kk) - ht*W*F_k ) ) )); 
end
toc

U = zeros(Ntot, Nt);
U(dof, :) = U_exp(:, :);

% recovering the physical variable Y = (Y1, Y1)^T
Y1 = 1/sqrt(2)*(U(NNB(2, :), :) - U(NNB(1, :), :));
Y2 = 1/sqrt(2)*(U(NNB(2, :), :) + U(NNB(1, :), :));

Y2ell = zeros(1, Nt);
for kk = 1:Nt
    Y2ell(1, kk) = Y2(Nx, kk);
end
figure()
plot(t, Y2ell(1, :));
xlabel('time');
title('y2 at x = ell');

tic
%% Plot of physical variable
%[taxis,xaxis] = meshgrid(time,x);
figure();

subplot(1, 2, 1);
s = surf(t,x, Y1);
s.EdgeColor = 'none';
xlabel('time'); ylabel('space'); zlabel('solution 1');
title('transport first equation y1');

subplot(1, 2, 2);
s = surf(t,x,Y2);
s.EdgeColor = 'none';
xlabel('time'); ylabel('space'); zlabel('solution 2');
title('transport second equation y2');
toc


%% Plot of diagonal variable
figure()

subplot(1, 2, 1);
s = surf(t,x,U(NNB(1, :), :));
s.EdgeColor = 'none';
xlabel('time'); ylabel('space'); zlabel('solution 1');
title('transport first equation r-');


subplot(1, 2, 2);
s = surf(t,x,U(NNB(2, :), :));
s.EdgeColor = 'none';
xlabel('time'); ylabel('space'); zlabel('solution 2');
title('transport second equation r+');
toc