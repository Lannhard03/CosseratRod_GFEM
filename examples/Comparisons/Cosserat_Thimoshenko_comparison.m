function [F, F_l] = test_network_local(m, angle)
F = FiberNetwork;

F.m = m;
F.G = 1/(F.m - 1) * ones(F.m - 1, 1); %Uniform grid


F.edges = zeros(F.m - 1, 2);
F.edges(:, 1) = 1:(F.m - 1);
F.edges(:, 2) = 2:F.m;

F.A = repmat([1963, 755, 755], F.m - 1, 1);
F.K = repmat([0.94, 1.23, 1.23], F.m - 1, 1);

F.fixed_nodes = [1, F.m];
F.connected_nodes = [];

F = F.apply_connections();

F.phi = zeros(F.m, 3);
F.R = zeros(F.m, 4);

F.phi(:, 1) = linspace(0, 1, F.m).';
%F.R(:, 1) = 1;

angle = linspace(0, angle, F.m)/2;
F.R(:, 1) = cos(angle);
F.R(:, 4) = sin(angle);

%Hard code R
%F.R(1, :) = [1, 0, 0, 0];
%F.R(F.m, :) = [1, 0, 0, 0];
F.U_pre = zeros(F.m - 1, 3); %One strain for each edge.
F.U_pre(:, 1) = 1;

F.R_pre = zeros(size(F.edges, 1), 3);

F = F.update_network();

F_l = normal2linearNetwork(F);
end

%Define solver parameters
prm = solver_param;
prm.backtracking_max = 100;
prm.beta = 0.9;
prm.alpha = 0.6;
prm.sigma = 1e-12;
prm.iter_max = 1000;
prm.allowed_failed_backtracks = 5;


m = 64;

num = 5;
angles = pi * 2.^(-(1:num));

dists = zeros(1, num);

for i = 1:num
    angle = angles(i);
    [F, F_l] = test_network_local(64, angle);

    F_l = thimoshenko_solver(F_l, prm);
    F   = exact_cosserat_solver(F, prm);
        
    phi_l = F_l.start + F_l.u;

    diff = vecnorm(F.phi - phi_l, 2, 2);
    dists(i) = max(diff);
end

%draw_network(F);
%draw_linearnetwork(F_l);
