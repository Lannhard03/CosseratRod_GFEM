function F = test_network_local(m, angle)
    F = LinearFiberNetwork;
    
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
    
    F.u = zeros(F.m, 3);
    
    F.start = zeros(F.m, 3);
    F.start(:, 1) = linspace(0, 1, F.m).';
    
    R = zeros(F.m, 4);
    angles = linspace(0, angle, F.m)/2;
    R(:, 1) = cos(angles);
    R(:, 4) = sin(angles);
    F.r = linearize_rotations(R);
    
    %Hard code R
    F.n = zeros(F.m - 1, 3); %One strain for each edge.
    F.n(:, 1) = 1;
end

function Rs = estimate_strain(F)
    Rs = zeros(size(F.edges, 1), 3);
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
    
        R1 = F.R(left, :);
        R2 = F.R(right, :);
    
        h = F.G(i);
    
        %Strain
        logR = quat_log(R1, R2);
        Rs(i, :) = logR(1) * logR(2:4) / h;
    end
end

%Cool values:
%1. m = 64, R_pre = [0.1, 0, 0], alpha = 0.6 or 1.0 (convergence of 0.6 but
%not for 1.0).

%Define solver parameters
prm = solver_param;
prm.backtracking_max = 100;
prm.beta = 0.9;
prm.alpha = 0.6;
prm.sigma = 1e-12;
prm.iter_max = 1000;
prm.allowed_failed_backtracks = 5;


m = 64;
F = test_network_local(64, pi/4);

draw_linearnetwork(F);

F = thimoshenko_solver(F, prm);

draw_linearnetwork(F);
