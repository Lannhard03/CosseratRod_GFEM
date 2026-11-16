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
F.r = zeros(F.m, 3);

%F.R(:, 1) = 1;

R = zeros(F.m, 4);
a = 0; 
b = pi/2; 
angle = linspace(a, b, 32)/2;
R(:, 1) = cos(angle);
R(:, 4) = sin(angle);
R(F.m, :) = [cos(pi/4), 0, 0, sin(pi/4)];
R(1, :) = [cos(pi/4), 0, sin(pi/4), 0];
F.r = linearize_rotations(R);


F.start = zeros(F.m, 3);
F.start(:, 1) = linspace(0, 0.5, 32).';

F.n = zeros(F.m - 1, 3); 
F.n(:, 1) = 1;

%Define solver parameters
prm = solver_param;
prm.backtracking_max = 1000;
prm.beta = 0.9;
prm.alpha = 2;
prm.sigma = 1e-4;
prm.iter_max = 1000;
prm.allowed_failed_backtracks = 5;

draw_linearnetwork(F);

F = thimoshenko_solver(F, prm);

draw_linearnetwork(F);

%Compare with reference solution
load("data/32grid_RefSol.mat");
%draw_rod(phi_ref, R_ref);

%phi = F.phi;

%rod1_len = approx_rod_len(F.phi(1:32, :))
%max_dist_to_ref = max(vecnorm(phi - phi_ref, 2, 2))

