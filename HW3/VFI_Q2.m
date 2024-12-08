beta = 0.95;
alpha = 0.4;
delta = 0.06;
u = @(c) log(c);
k_min = 0;
k_max = 10;
Nk = 100;
k_grid = linspace(k_min, k_max, Nk)';

Nz = 5;
mu = 0;
rho = 0.9;
sigma = 0.1;
m = 1;
[Z,Zprob] = tauchen(Nz,mu,rho,sigma,m);

V = zeros(Nk, Nz);
policy_k = zeros(Nk, Nz);

max_iter = 1000; 
tol = 1e-6;    
for iter = 1:max_iter
    V_new = zeros(Nk, Nz);
    for i = 1:Nk
        for j = 1:Nz
            z0 = Z(j); 
            k0 = k_grid(i);

            c = exp(z0) * k0^alpha + (1 - delta) * k0 - k_grid;
            U = u(c);
            U(c <= 0) = -inf; 

            EV = V * Zprob(j,:)'; 

            total_value = U + beta * EV;

            [V_new(i,j), policy_index] = max(total_value);
            policy_k(i,j) = k_grid(policy_index);
        end
    end
    
    if max(abs(V_new(:) - V(:))) < tol
        disp(['Converged (number of iterations: ', num2str(iter), ')']);
        break;
    end
    V = V_new;
end

figure;
hold on;
for i_z = 1:Nz
    plot(k_grid, V(:, i_z), 'DisplayName', ['z = ', num2str(Z(i_z))]);
end
xlabel('Capital k');
ylabel('Value function V(k, z)');
title('Value functions (for each z_t)');
legend show;
grid on;

figure;
hold on;
for i_z = 1:Nz
    plot(k_grid, policy_k(:, i_z), 'DisplayName', ['z = ', num2str(Z(i_z))]);
end
xlabel('Capital k');
ylabel('Optimal next capital k''');
title('Policy functions (for each z_t)');
legend show;
grid on;

