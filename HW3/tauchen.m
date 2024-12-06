function [Z,Zprob] = tauchen(N,mu,rho,sigma,m)

N = 5; % Number of Grid Points, the number of potential realizations of z.
mu = 0; % Mean
rho = 0.9; % AR(1) Coefficient
sigma = 0.1; % Standard Deviation
m = 3; % Number of Standard Deviations
    
Z = zeros(N,1); % Grid
Zprob = zeros(N,N); % Transition Matrix
c = (1-rho)*mu; % Constant

% Define Grids 
zmax  = m*sqrt(sigma^2/(1-rho^2));
zmin  = -zmax;
w = (zmax-zmin)/(N-1);
Z = linspace(zmin,zmax,N)';

% Stationary value, mu
Z = Z + mu;

% Create Transition Matrix
for j = 1:N
    for k = 1:N
        if k == 1
            Zprob(j,k) = normcdf((Z(1)-c-rho*Z(j)+w/2)/sigma);
        elseif k == N
            Zprob(j,k) = 1 - normcdf((Z(N)-c-rho*Z(j)-w/2)/sigma);
        else
            Zprob(j,k) = normcdf((Z(k)-c-rho*Z(j)+w/2)/sigma) - ...
                         normcdf((Z(k)-c-rho*Z(j)-w/2)/sigma);
        end
    end
end