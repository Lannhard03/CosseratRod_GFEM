F = FiberNetwork;

%Start of rod 1 and end of rod 2 are fixed
F.fixed_nodes = [1, 64]; 
%Nodes 32 and 33 are connected: their relative positions and rotations should be preserved.
%For more connections add more rows.
F.connected_nodes = [32, 33]; 

%Define initial state of network
rod1 = Rod;
rod2 = Rod;

rod1.m = 32;
rod1.G = linspace(0,1, rod1.m);
rod1.A = [1963, 755, 755];
rod1.K = [0.94, 1.23, 1.23];

rod2.m = 32;
rod2.G = linspace(0,1, rod2.m);
rod2.A = [1963, 755, 755];
rod2.K = [0.94, 1.23, 1.23];

rod1.phi = zeros(rod1.m, 3);
rod1.R = zeros(rod1.m, 4);

rod1.phi(:, 1) = linspace(0, 1, 32).';
rod1.R(:, 1) = 1;

rod2.phi = zeros(rod2.m, 3);
rod2.R = zeros(rod2.m, 4);

rod2.phi(:, 1) = 1;
%Add "bending" here!
rod2.phi(:, 1) = linspace(1, 0.5, 32).';
rod2.phi(:, 2) = linspace(0, 1, 32).';

rod2.R(:, 1) = 1/sqrt(2);
rod2.R(:, 4) = 1/sqrt(2); %Rotation 90 around z-axis

rod1.Us_pre = zeros(rod1.m, 3);
rod1.Us_pre(:, 1) = 1;

rod2.Us_pre = zeros(rod2.m, 3);
rod2.Us_pre(:, 1) = 1;

rod1.Rs_pre = zeros(rod1.m, 3);
rod2.Rs_pre = zeros(rod2.m, 3);

F.rods = [rod1, rod2];

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

F = BasicSolver(F, prm);

%Draw "solved" rod
draw_network(F)