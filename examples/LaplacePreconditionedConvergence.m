function F = test_network(m)
    F = FiberNetwork;
    
    F.m = m;
    F.G = 1/(F.m - 1) * ones(F.m - 1, 1); %Uniform grid
    
    F.edges = zeros(F.m - 1, 2);
    F.edges(:, 1) = 1:(F.m - 1);
    F.edges(:, 2) = 2:F.m;
    
    F.fixed_nodes = [1, F.m];
    F.connected_nodes = [];
    
    F = F.apply_connections();
    
    F.phi = zeros(F.m, 3);
    F.R = zeros(F.m, 4);
    
    %Phi goes forward for visualization purposes, we only solve over R!!
    F.phi(:, 1) = linspace(0, 1, F.m).';
    %F.R(:, 1) = 1;
    
    a = 0; 
    b = pi/2; 
    angle = linspace(a, b, F.m)/2;
    F.R(:, 1) = cos(angle);
    F.R(:, 3) = sin(angle);
    
    %Hard code R
    F.R(1, :) = [1, 0, 0, 0];
    F.R(F.m, :) = [cos(pi/4), 0, 0, sin(pi/4)];
end


%Define solver parameters
prm = solver_param;
prm.backtracking_max = 100;
prm.beta = 0.9;
%2 / (1 + (1 - dist/8)) = 1.1612 is optimal for Richardson iteration,
% but 2 works better here? Maybe bad estimate of eigenvals?
prm.alpha = 1.1612; 
prm.sigma = 1e-4;
prm.iter_max = 1000;
prm.allowed_failed_backtracks = 5;


%Exact energy is 1/2 * (dist(u_0, u_1) * 1/sqrt(2))^2 = pi^2 / 16
exact_energy = pi^2 / 16;

%Theoretical convergence rate is 1/8 * dist(u_0, u_1) = 
% = sqrt(2)pi/16 or roughly 0.2777
conv_rate_theo = sqrt(2)*pi/16;

energy_hist = zeros(4, prm.iter_max);
m = [16, 32, 64, 128, 256];
for i = 1:5
    F = test_network(m(i));
    [F, energy_hist(i, :)] = laplace_solver(F, prm);
end

energies = energy_hist(5, :);
num = sum(energies ~= 0);
diff = abs(energies(1:num) - exact_energy);
convergence_rate = diff(2:num) ./ diff(1:(num-1));

hold on
plot(1:num, diff)
plot(1:num, conv_rate_theo.^(1:num))
yscale log