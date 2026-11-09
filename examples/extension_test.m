F = FiberNetwork;

F.m = 32;
F.G = 1/(F.m - 1) * ones(F.m - 1, 1); %Uniform grid

F.edges = zeros(F.m - 1, 2);
F.edges(:, 1) = 1:31;
F.edges(:, 2) = 2:32;

F.fixed_nodes = [1, 32];
F.connected_nodes = {};

F = F.apply_connections();

F.A = repmat([1963, 755, 755], F.m - 1, 1);
F.K = repmat([0.94, 1.23, 1.23], F.m - 1, 1);

F.phi = zeros(F.m, 3);
F.R = zeros(F.m, 4);

F.phi(:, 1) = linspace(0, 1, 32).';
F.phi(32, 1) = 0.5;
F.R(:, 1) = 1;

a = 0; 
b = pi/2; 
angle = linspace(a, b, 32)/2;
F.R(:, 1) = 1;

F.U_pre = zeros(F.m - 1, 3); %One strain for each edge.
F.U_pre(:, 1) = 1;

F.R_pre = zeros(F.m - 1, 3);

draw_network(F);

%Define solver parameters
prm = solver_param;
prm.backtracking_max = 1000;
prm.beta = 0.9;
prm.alpha = 2;
prm.sigma = 1e-4;
prm.iter_max = 1000;
prm.allowed_failed_backtracks = 5;

F =  cosserat_solver(F, prm);

draw_network(F);


