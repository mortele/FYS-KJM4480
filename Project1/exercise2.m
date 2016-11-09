function [ ] = exercise2()
close all;
clear variables;
clc;
format;


%% Parameters
N = 2;      % Number of particles.


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
% the spin-orbital index, i.
e = @(i) 1*(i>=1 & i<=2) + 2*(i>=3 & i<=6) + 3*(i>=7 & i<=12);

% Spin values of the spin-orbitals labeled by the index i.
spin = @(i) -1 + (mod(i,2)~=0)*2;

% Map the spin-orbitals onto the corresponding orbitals, using the
% numbering 1, 2  --> (0,0)  = 1
%           3, 4  --> (0,-1) = 2 
%           5, 6  --> (0,+1) = 3 
%           7, 8  --> (0,-2) = 4 
%           9, 10 --> (1,0)  = 5 
%           11,12 --> (0,+2) = 6.
mapState = @(i) (i+mod(i,2))/2;

% w-integral value of the mapped states.
wm = @(p,q,r,s) w(mapState(p), mapState(q), mapState(r), mapState(s));

% Anti-symmetrized w.
function [w] = wAS(p,q,r,s)
    wpqrs = delta(spin(p),spin(r)) * delta(spin(q),spin(s)) * wm(p,q,r,s);
    wpqsr = delta(spin(p),spin(s)) * delta(spin(q),spin(r)) * wm(p,q,s,r);
    w = wpqrs - wpqsr;
end


%% Compute the expectation value of the Hamiltonian, Eref
E_ref   = 0; % Expectation value of the total Hamiltonian.

% One-body operator contribution to the Hamiltonian, H0.
for i=1:N
    E_ref = E_ref + e(i);
end

% Two-body operator contribution to the Hamiltoninan, W.
for i=1:N
    for j=1:N
        E_ref = E_ref + 0.5 * wAS(i,j,i,j);
    end
end


%% Set up the CIS matrix
H = zeros(3,3);

% Elements of the CIS matrix.
H(1,1) = E_ref;
for i=1:N
    H(1,2) = H(1,2) + wAS(1,i,9,i);
    H(1,3) = H(1,3) + wAS(2,i,10,i);
end
H(2,3) = 3 * wAS(2,9,10,1) / 4;

H(2,2) = E_ref - e(1) + e(9)  + 3 * wAS(1, 9, 9,1) / 4;
H(3,3) = E_ref - e(2) + e(10) + 3 * wAS(2,10,10,2) / 4;
for i=1:N
    H(2,2) = H(2,2) - wAS(1,i,1,i) + wAS(9, i, 9,i);
    H(3,3) = H(3,3) - wAS(2,i,2,i) + wAS(10,i,10,i);
end

% Copy the upper diagonal onto the lower diagonal of CIS, and diagonalize
% the matrix.
H = H + triu(H,1)';
E_CIS = eig(H);

fprintf('E_ref = %f, E_CIS = %f\n\n',E_ref, E_CIS(1));

assignin('base', 'A', A);
assignin('base', 'E_ref', E_ref);
assignin('base', 'H', H);
assignin('base', 'E_CIS', E_CIS);
exit;
end