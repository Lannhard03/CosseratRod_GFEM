F = FiberNetwork;

F.m = 32;
F.G = 1/(F.m - 1) * ones(F.m - 1, 1); %Uniform grid

F.edges = zeros(F.m - 1, 2);
F.edges(:, 1) = 1:31;
F.edges(:, 2) = 2:32;

F.fixed_nodes = [1, 32];
F.connected_nodes = [];

F = F.apply_connections();

F.A = repmat([1963, 755, 755], F.m - 1, 1);
F.K = repmat([0.94, 1.23, 1.23], F.m - 1, 1);

F.phi = zeros(F.m, 3);
F.R = zeros(F.m, 4);

F.phi(:, 1) = linspace(0, 0.5, 32).';
F.R(:, 1) = 1;

a = 0; 
b = pi/2; 
angle = linspace(a, b, 32)/2;
F.R(:, 1) = cos(angle);
F.R(:, 4) = sin(angle);
F.R(F.m, :) = [cos(pi/4), 0, 0, sin(pi/4)];
F.R(1, :) = [cos(pi/4), 0, sin(pi/4), 0];

F.Us_pre = zeros(F.m - 1, 3); %One strain for each edge.
F.Us_pre(:, 1) = 1;

F.Rs_pre = zeros(F.m - 1, 3);

draw_network(F);


%Define solver parameters
prm = solver_param;
prm.backtracking_max = 1000;
prm.beta = 0.9;
prm.alpha = 0.1;
prm.sigma = 1e-4;
prm.iter_max = 1000;
prm.allowed_failed_backtracks = 5;

F = BasicSolver(F, prm);

draw_network(F);

%Compare with reference solution
load("data/32grid_RefSol.mat");
draw_rod(phi_ref, R_ref);

phi = F.phi;

rod1_len = approx_rod_len(F.phi(1:32, :))
max_dist_to_ref = max(vecnorm(phi - phi_ref, 2, 2))

