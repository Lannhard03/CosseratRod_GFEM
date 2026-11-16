F = LinearFiberNetwork;

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

F.u = zeros(F.m, 3);
F.start = zeros(F.m, 3);
F.r = zeros(F.m, 3);

F.start(:, 1) = linspace(0, 1, 32).';

F.u(32, 1) = -0.5;

F.n = zeros(F.m - 1, 3); %One strain for each edge.
F.n(:, 1) = 1;

draw_linearnetwork(F);

%Define solver parameters
prm = solver_param;
prm.backtracking_max = 1000;
prm.beta = 0.9;
prm.alpha = 2;
prm.sigma = 1e-4;
prm.iter_max = 1000;
prm.allowed_failed_backtracks = 5;

F = thimoshenko_solver(F, prm);

draw_linearnetwork(F);


