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


%% Functions
delta   = @(a,b)            a == b;
H0      = @(p1,q1,p2,q2)    2*(p1-1)*delta(p1,p2)*delta(q1,q2) + ...
                            2*(q1-1)*delta(p1,p2)*delta(q1,q2);
V       = @(p1,q1,p2,q2,g)  -0.5*g*(   delta(p1,p2) + delta(p1,q2) ...
                                     + delta(q1,p2) + delta(q1,q2)    );
H       = @(p1,q1,p2,q2,g) H0(p1,q1,p2,q2) + V(p1,q1,p2,q2,g);
map     = @(I) [1,2]*(I==1) + [1,3]*(I==2) + [1,4]*(I==3) + ...
               [2,3]*(I==4) + [2,4]*(I==5) + [3,4]*(I==6);
           
%% Set up and diagonalize FCI matrix
FCI = zeros(basisSize, basisSize);
G   = linspace(-1,1,n); 
E   = zeros(n,basisSize);
f   = zeros(n,1);

for k=1:n
    g = G(k);
    for i=1:basisSize
        for j=1:basisSize
            p1q1 = map(i);
            p2q2 = map(j);
            
            p1   = p1q1(1);
            q1   = p1q1(2);
            
            p2   = p2q2(1);
            q2   = p2q2(2);
            
            FCI(i,j) = H(p1,q1,p2,q2,g);
        end
    end
    if ~ishermitian(FCI) 
        error = 'FCI matrix is not hermitian.'
        FCI
        break
    end
    [V, e]      = eig(FCI);
    e           = diag(e);         
    E(k,:)      = e;
    GS          = zeros(basisSize,1); GS(1) = 1;
    f(k)        = V(1,1)^2;
end


%% Plot the energies and the prop. of unperturbed GS
% Energy plot.
figure(1);
legendA = cell(basisSize,1);
for k=1:basisSize
    plot(G, E(:,k));
    legendA{k} = sprintf('E%d', k);
    hold on;
end
legend(legendA)
xlabel('g');
ylabel('E');

% Probability of finding the perturbed system in the unperturbed ground
% state, |12> (index 1 under the mapping used here).
figure(2);
plot(G, f);
xlabel('g');
ylabel('f_0');

    
    
    
    
    
    
    
    
