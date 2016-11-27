function [ ] = exercise3( )
%% Clean up
close all;
clear all;
clc;
format long;


%% Numerics and constants
eps             = 1e-15;        % Convergence criterion.
sigma           = -0.5;         % Convergence helper factor.
maxIterations   = 100;          % Maximum number of iterations on G.


%% Functions
% Fock operator values, f_p
f   = @(p,g) (p-1) + (p==1 || p==2)*(-g/2);

% Map onto the 2 x 2 - matrix of t amplitudes.
map = @(p)   (p==1 || p==3)*1 + (p==2 || p==4)*2;

% The function we are iterating on, G(t).
function [t_] = G(t, i, a, sigma, g)
    denominator = 1/(2 * (f(i,g) - f(a,g)) + sigma);
    numerator   = sigma * t(map(i), map(a)) ...
                    - g/2 * (1 + t(map(i), map(3)) + t(map(i), map(4)) ...
                               + t(map(1), map(a)) + t(map(2), map(a)) ...
                               + t(map(1), map(3)) * t(map(2), map(4)) ...
                               + t(map(1), map(4)) * t(map(2), map(3)));
    t_ = numerator * denominator;
end

% The CCD energy.
function [E_CCD] = E(t, g) 
    singleStateEnergy1 = 0;
    singleStateEnergy2 = 1;
    E_CCD = 2*singleStateEnergy1 + 2*singleStateEnergy2 -g ...
                -g/2 * sum(sum(t));
end

%% Perform the CCD iterations.
N  = 100;
g_ = linspace(-1,1,N);
e_ = zeros(N,1);

% Iterate on G.
for i=1:N
    % Update value of g.
    g = g_(i);
    
    % Initial guess for the t ampltitudes is just a zero matrix.
    t  = zeros(2,2);
    t_ = zeros(2,2);
    
    for iteration=1:maxIterations
        % Perform the iteration on the function G.
        t_(1,1) = G(t, 1, 3, sigma, g);
        t_(1,2) = G(t, 1, 4, sigma, g);
        t_(2,1) = G(t, 2, 3, sigma, g);
        t_(2,2) = G(t, 2, 4, sigma, g);
        
        % Check if we have reached self-consistency.
        if max(max(abs(t-t_))) < eps 
            break;
        elseif iteration == maxIterations
            error = 'Self-consistency not reached for'
            g
        end
        
        % Update the t matrix with the new values.
        t = t_;
    end
    
    e_(i,1) = E(t_, g);
end

save('mongo', 'e_', 'g_');
%% Plot the resulting energy
figure(1);
plot(g_, e_, 'r-');
xlabel('$g$', 'FontSize', 16, 'interpreter', 'latex');
ylabel('G-S. energy', 'FontSize', 16, 'interpreter', 'latex');
leg = legend('CCD');
set(leg, 'FontSize', 16, 'interpreter', 'latex');

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height]; 

end