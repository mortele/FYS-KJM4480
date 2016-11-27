%% Clean up
close all;
clear variables;
clc;
format;

%% Constants / numerics
M           = 4;
N           = 4;
basisSize   = 6;
n           = 100;
plotting    = true;

%% Functions
% Simple delta function.
delta   = @(a,b) a == b;

% One body operator expectation value, <p1q1|H0|p2q2>.
H0      = @(p1,q1,p2,q2)    2*(p1-1)*delta(p1,p2)*delta(q1,q2) + ...
                            2*(q1-1)*delta(p1,p2)*delta(q1,q2);

% Two body operator expectation value, <p1q1|V|p2q2>.
V       = @(p1,q1,p2,q2,g)  -0.5*g*(   delta(p1,p2) + delta(p1,q2) ...
                                     + delta(q1,p2) + delta(q1,q2)    );
                                 
% Full Hamiltonian expectation value, <p1q1|H|p2q2>.
H       = @(p1,q1,p2,q2,g) H0(p1,q1,p2,q2) + V(p1,q1,p2,q2,g);

% Map two indices (a state |ij>) onto a single index.
map     = @(I) [1,2]*(I==1) + [1,3]*(I==2) + [1,4]*(I==3) + ...
               [2,4]*(I==4) + [2,3]*(I==5) + [3,4]*(I==6);

% The constant prefactor of the operator R (1/(e0-ej)), in terms of a
% single index (as mapped above).
gamma   = @(I) 1./(H0(1,2,1,2) - (I==2) * H0(1,3,1,3) ...
                               - (I==3) * H0(1,4,1,4) ...
                               - (I==4) * H0(2,4,2,4) ...
                               - (I==5) * H0(2,3,2,3) ...
                               - (I==6) * H0(3,4,3,4));

%% Set up and diagonalize FCI matrix
FCI     = zeros(basisSize,   basisSize);
CID     = zeros(basisSize-1, basisSize-1);
G       = linspace(-1,1,n); 
E_FCI   = zeros(n,basisSize);
E_CID   = zeros(n,basisSize-1);
E_RS1   = zeros(n,1);
E_RS2   = zeros(n,1);
E_RS3   = zeros(n,1);
f_FCI   = zeros(n,1);
f_CID   = zeros(n,1);

for k=1:n
    g = G(k);
    for i=1:basisSize
        for j=1:basisSize
            p1q1     = map(i); p1 = p1q1(1); q1 = p1q1(2);
            p2q2     = map(j); p2 = p2q2(1); q2 = p2q2(2);
            FCI(i,j) = H(p1,q1,p2,q2,g);
        end
    end
    
    % Set up the CID matrix.
    CID = FCI(1:end-1, 1:end-1);
    
    % Check that the FCI and the CID matrices are hermitian (as they should
    % be).
    if ~ishermitian(FCI) || ~ishermitian(CID)
        error = 'FCI or CID matrices are not hermitian.'
        FCI
        CID
        break
    end
    
    % Diagonalize the FCI and the CID matrices, extract the ground state
    % energy and the ground state probability, f.
    [V_FCI, e_FCI]  = eig(FCI);
    [V_CID, e_CID]  = eig(CID);
    e_FCI           = diag(e_FCI);
    e_CID           = diag(e_CID);
    E_FCI(k,:)      = e_FCI;
    E_CID(k,:)      = e_CID;
    GS_FCI          = [1 0 0 0 0 0];
    GS_CID          = [1 0 0 0 0];
    f_FCI(k)        = V_FCI(1,1)^2;
    f_CID(k)        = V_CID(1,1)^2;
    
    % Set up the Rayleigh-Schrödinger pertubation theory to third order.
    E0  = H0(1,2,1,2);
    E1  = -g;
    E2  = -g/2 * (gamma(2)*FCI(1,2) + gamma(3)*FCI(1,3) + ...
                  gamma(4)*FCI(1,4) + gamma(5)*FCI(1,5));
    E3A = -g^2/4 * ((2*gamma(2)^2 + gamma(2)*gamma(3) + gamma(2)*gamma(4))*FCI(1,2) + ...
                    (2*gamma(3)^2 + gamma(2)*gamma(3) + gamma(3)*gamma(5))*FCI(1,3) + ...
                    (2*gamma(4)^2 + gamma(2)*gamma(3) + gamma(4)*gamma(5))*FCI(1,4) + ...
                    (2*gamma(5)^2 + gamma(3)*gamma(5) + gamma(4)*gamma(5))*FCI(1,5) + ...
                    (gamma(2)*gamma(6) + gamma(3)*gamma(6) + gamma(4)*gamma(6) + gamma(5)*gamma(6))*FCI(1,6)) + ...
          -g^2/2 * ( gamma(2)^2*FCI(2) + gamma(3)^2*FCI(3) + gamma(4)^2*FCI(4) + gamma(5)^2*FCI(5));          
    E3B = E1 * (-g)*(gamma(2)^2*FCI(1,2) + gamma(3)^2*FCI(1,3) + ...
                     gamma(4)^2*FCI(1,4) + gamma(5)^2*FCI(1,5));
    E3  = E3A + E3B;
    
    e_RS1      = E0 + E1;
    e_RS2      = e_RS1 + E2;
    e_RS3      = e_RS2 + E3;
    E_RS1(k,1) = e_RS1;
    E_RS2(k,1) = e_RS2;
    E_RS3(k,1) = e_RS3;
end

%% Plot the energies and the prop. of unperturbed GS
% Energy plot (FCI).
if plotting 
    figure(1);
    legendA = cell(basisSize,1);
    for k=1:basisSize
        if k==1
            plot(G, E_FCI(:,k), 'b-', 'LineWidth', 2);
        else 
            if mod(k,2)==0
                plot(G, E_FCI(:,k), '--');
            else 
                plot(G, E_FCI(:,k), '-');
            end
        end
        legendA{k} = sprintf('E%d', k);
        hold on;
    end
    leg = legend(legendA);
    set(leg, 'FontSize', 16, 'Interpreter', 'latex');
    xlabel('g', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('Energy(g)', 'FontSize', 16, 'Interpreter', 'latex');
    title('Full configuration interaction');

    % Energy plot (CID).
    figure(2);
    legendA = cell(basisSize-1,1);
    for k=1:basisSize-1
        if k==1
            plot(G, E_CID(:,k), 'b-', 'LineWidth', 2);
        else 
            if mod(k,2)==0
                plot(G, E_CID(:,k), '--');
            else 
                plot(G, E_CID(:,k), '-');
            end
        end
        legendA{k} = sprintf('E%d', k);
        hold on;
    end
    leg = legend(legendA);
    set(leg, 'FontSize', 16, 'Interpreter', 'latex');
    xlabel('$g$', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('Energy($g$)', 'FontSize', 16, 'Interpreter', 'latex');
    %title('Configuration interaction doubles');
    
    % Energy plot (RS3).
    figure(3);
    plot(G, E_RS3);
    xlabel('g');
    ylabel('E0');
    title('Rayleigh-Schrödinger 3');


    % Probability of finding the perturbed system in the unperturbed ground
    % state, |12> (index 1 under the mapping used here).
    figure(4);
    hold('on')
    %plot(G, f_FCI, '-');
    plot(G, f_CID, '--');
    xlabel('$g$', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('Probability of GS, $f_0$', 'FontSize', 16, 'Interpreter', 'latex');
    leg = legend('CID');%,'CID');
    set(leg, 'FontSize', 16, 'Interpreter', 'latex');
    
    figure(5);
    hold('on');
    plot(G, E_RS1, 'b--');
    plot(G, E_RS2, 'r--');
    plot(G, E_RS3, 'k-');
    xlabel('$g$', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('G-S. Energy($g$)', 'FontSize', 16, 'Interpreter', 'latex');
    leg = legend('RSPT1','RSPT2','RSPT3');
    set(leg, 'FontSize', 16, 'Interpreter', 'latex');
    
    figure(6);
    hold('on');
    plot(G, E_FCI(:,1), 'b-');
    plot(G, E_CID(:,1), 'r-');
    plot(G, E_RS3, 'k-');
    xlabel('$g$', 'FontSize', 16, 'Interpreter', 'latex');
    ylabel('G-S. Energy($g$)', 'FontSize', 16, 'Interpreter', 'latex');
    leg = legend('FCI','CID','RSPT3');
    set(leg, 'FontSize', 16, 'Interpreter', 'latex');
    
end


%%
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];    
    
    
