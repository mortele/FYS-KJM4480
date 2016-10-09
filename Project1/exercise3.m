function [ ] = exercise3()
close all;
clear variables;
clc;
format short;


%% Parameters
N               = 2;  % Number of particles.
L               = 12; % Number of basis functions. 
maxIterations   = 20; % Maximum number of SCT iterations.


%% Load integrals from file
fileName    = 'coulomb.dat';
inFile      = fopen(fileName, 'rt');
integrals   = textscan(inFile, '%f %f %f %f %f');
A           = cell2mat(integrals);
fclose(inFile);


%% Functions
% Function for extracting the value of the w-integral from the
% 'coulomb.dat' file.
w = @(p,q,r,s) A(A(:,1)==p & A(:,2)==q & A(:,3)==r & A(:,4)==s, 5);

% Kroenecker delta.
delta = @(i,j) i==j;

% Values of the single-particle part of the Hamiltonian, as a function of
% the spatial index, i.
e = @(i) 1*(i==1) + 2*(i>=2 & i<=3) + 3*(i>=4 & i<=6);

% Computes [qr|ps].
function [W] = qrps(q,r,p,s)
    qrps = w(q,r,p,s);
    qrsp = w(q,r,s,p);
    W = 0;
    if isempty(qrps)==false
        W = W + qrps;
    end
    if isempty(qrsp)==false
        W = W - 0.5 * qrsp;
    end
end

% Density matrix.
function [D] = densityMatrix(U) 
    D = zeros(L/2, L/2);
    for r=1:L/2
        for s=1:L/2
            D(s,r) = 0;
            for j=1:N/2
                D(s,r) = D(s,r) + 2.0 * U(s,j) * conj(U(r,j));
            end
        end
    end
end
                    
% Fock matrix.
function [F] = FockMatrix(D)
    F = diag([1 2 2 3 3 3]');
    for q=1:L/2
        for p=1:L/2
            for r=1:L/2
                for s=1:L/2
                    F(q,p) = F(q,p) + D(s,r) * qrps(q,r,p,s);
                end
            end
        end
    end
end

% Computes the restricted Hartree-Fock energy, given a diagonalized Fock
% matrix.
function [E_RHF] = restrictedHartreeFockEnergy(epsilon, U, D)
    E_RHF = 0;
    W = 0;
    for i=1:N/2
        E_RHF = E_RHF + 2 * epsilon(i,i);
        for q=1:L/2
            for p=1:L/2
                srSum = 0;
                for r=1:L/2
                    for s=1:L/2
                        srSum = srSum + D(s,r)*qrps(q,r,p,s);
                    end
                end
                W = W + conj(U(q,i)) * srSum * U(p,i);
            end
        end
    end
    E_RHF = E_RHF - W;
end    

% Checks if a given matrix is Hermitian or not.
function [hermitian] = checkHermitian(A)
    if max(max(A-A'))>eps
        hermitian = false;
    else
        hermitian = true;
    end
end


%% Perform SCF iterations
% Initial guess for the unitary U matrix is I.
U = eye(L/2,L/2);

% Vector containing the progression of HF-energy during the SCF iterations.
E_SCF  = zeros(maxIterations,1);
deltak = zeros(maxIterations,1);

% Plot the progress live.
figure(1);
subplot(2,1,2);
h2 = semilogy(nan, nan, 'r-o');
xlabel('Iteration', 'interpreter', 'latex','FontSize', 16);
ylabel('$\delta_k$', 'interpreter', 'latex','FontSize', 16);
subplot(2,1,1);
h1 = plot(nan, nan, 'r-o');
ylabel('$E_{RHF}$', 'interpreter', 'latex','FontSize', 16);
oldEpsilon = 1;

for iteration=1:maxIterations
    % Compute D and F, and check that F is hermitian.
    D = densityMatrix(U);
    F = FockMatrix(D);
    H = checkHermitian(F);
    
    % Solve the linearized problem.
    [U, epsilon] = eig((F+F')/2);
    
    % Sort the eigenvalues and U.
    [e, I]  = sort(diag(epsilon), 'ascend');
    epsilon = epsilon(I,I);
    U       = U(:,I);
    
    % Store the current value of the energy.
    E_SCF(iteration)  = restrictedHartreeFockEnergy(epsilon, U, D);
    deltak(iteration) = max(abs(diag(epsilon)-diag(oldEpsilon)));
    
    % Update the plot.
    set(h1, 'XData', linspace(1,iteration,iteration));
    set(h1, 'YData', E_SCF(1:iteration));
    
    set(h2, 'XData', linspace(1,iteration,iteration));
    set(h2, 'YData', deltak(1:iteration));
    drawnow;
    
    % Print info to terminal.
    if H
        Hstr = 'Yes';
    else
        Hstr = 'No';
    end
    fprintf('Iteration: %2d  E = %-11.7f  F^H==F? %s  delta_k = %-11.7g\n',...
            iteration, ...
            E_SCF(iteration), ...
            Hstr,...
            deltak(iteration));
    oldEpsilon = epsilon;
    if abs(deltak(iteration)) < eps 
        break;
    end
end
end