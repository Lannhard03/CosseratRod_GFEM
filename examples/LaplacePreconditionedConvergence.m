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
prm.alpha = 2.0;
prm.sigma = 1e-4;
prm.iter_max = 1000;
prm.allowed_failed_backtracks = 5;

energy_hist = zeros(4, prm.iter_max);
m = [16, 32, 64, 128, 256];
for i = 1:5
    F = test_network(m(i));
    [F, energy_hist(i, :)] = LaplaceSolver(F, prm);
end


plot(1:prm.iter_max, energy_hist)
yscale log