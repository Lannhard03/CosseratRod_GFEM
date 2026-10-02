F = FiberNetwork;

F.m = 96;

num_edges = F.m - 3;
F.G = 3/(num_edges) * ones(num_edges, 1); %Uniform grid, in total length is 3.

F.edges = zeros(num_edges, 2);
F.edges(1:31, 1) = 1:31;
F.edges(1:31, 2) = 2:32;
F.edges(32:62, 1) = 33:63;
F.edges(32:62, 2) = 34:64;
F.edges(63:93, 1) = 65:95;
F.edges(63:93, 2) = 66:96;

F.fixed_nodes = [1, 64, 96];

F.connected_nodes = [32, 33, 65];
F = F.apply_connections();

F.A = repmat([1963, 755, 755], num_edges, 1);
F.K = repmat([0.94, 1.23, 1.23], num_edges, 1);

F.phi = zeros(F.m, 3);
F.R = zeros(F.m, 4);

F.phi(1:32, 1) = linspace(0, 1, 32).';
F.R(1:32, 1) = 1;

F.phi(33:64, 1) = linspace(1, 0.5, 32).';
F.phi(33:64, 2) = linspace(0, 1, 32).';

F.R(33:64, 1) = 1/sqrt(2);
F.R(33:64, 4) = 1/sqrt(2); %Rotation 90 around z-axis

F.phi(65:96, 1) = linspace(1, 0.5, 32).';
F.phi(65:96, 2) = linspace(0, -1, 32).';
F.phi(65:96, 3) = linspace(0, 0.5, 32).';

F.R(65:96, 1) = 1/sqrt(2);
F.R(65:96, 4) = -1/sqrt(2); %Rotation 90 around z-axis

F.Us_pre = zeros(num_edges, 3); %One strain for each edge.
F.Us_pre(:, 1) = 1;

F.Rs_pre = zeros(num_edges, 3);

draw_network(F);

init_log = quat_log(F.R(32, :), F.R(33, :));

%Define solver parameters
prm = solver_param;
prm.backtracking_max = 1000;
prm.beta = 0.9;
prm.alpha = 0.1;
prm.sigma = 1e-4;
prm.iter_max = 1000;
prm.allowed_failed_backtracks = 5;

F = BasicSolver(F, prm);


after_log = quat_log(F.R(32, :), F.R(33, :));

init_log
after_log

rod1_len = approx_rod_len(F.phi(1:32, :))
rod2_len = approx_rod_len(F.phi(33:64, :))
rod3_len = approx_rod_len(F.phi(65:96, :))

draw_network(F);