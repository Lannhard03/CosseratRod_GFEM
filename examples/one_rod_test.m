F = FiberNetwork;

%Start and end of rod are fixed
F.fixed_nodes = [1, 32]; 
%No connected nodes
F.connected_nodes = []; 

%Define initial configuration and pre-strain of material
rod1 = Rod;

rod1.m = 32;
rod1.G = linspace(0,1, rod1.m);
rod1.A = [1963, 755, 755];
rod1.K = [0.94, 1.23, 1.23];

rod1.phi = zeros(rod1.m, 3);
rod1.R = zeros(rod1.m, 4);

rod1.phi(:, 1) = linspace(0, 0.5, 32).';
rod1.R(:, 1) = 1;

a = 0; 
b = pi/2; 

angles = linspace(a/2, b/2, rod1.m);
rod1.R = zeros(rod1.m, 4);
rod1.R(:, 1) = cos(angles);
rod1.R(:, 4) = sin(angles);

rod1.R(1, :) = [cos(pi/4), 0, sin(pi/4), 0];

rod1.Us_pre = zeros(rod1.m, 3);
rod1.Us_pre(:, 1) = 1;

rod1.Rs_pre = zeros(rod1.m, 3);

F.rods = [rod1];

%Draw initial rod
draw_network(F);

%Define solver parameters
prm = solver_param;
prm.backtracking_max = 1000;
prm.beta = 0.9;
prm.alpha = 0.1;
prm.sigma = 1e-4;
prm.iter_max = 1000;
prm.allowed_failed_backtracks = 5;

%Run solver
F = BasicSolver(F, prm);

%Draw "solved" rod
draw_network(F);

%Compare with reference solution
load("data/32grid_RefSol.mat");
draw_rod(phi_ref, R_ref);

phi = F.rods(1).phi;

diff = max(vecnorm(phi - phi_ref, 2, 2))