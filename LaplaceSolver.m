%Solver for the "harmonic energy" problem:
%For u:[0,1] \to SO(3), minimize:
%\int_[0,1] <u', u'> dt.
%Exact solution is just a geodesic, so using GFEM is unnecessary,
%but good test problem.

function F = laplace_solver(F, params)

    J0 = J(F);    
    
    %h = F.G(1); 
    %D_h = diag(ones(F.m, 1)) - diag(ones(F.m - 1, 1), -1);
    %D_h = 1/h*D_h(:, 1:end-1);

    disc_laplacian = laplacian(F);
    
    %===== Outer loop =====
    for p = 1:params.iter_max
        fprintf("Iteration %d\n", p);
    
        gradR = grad_J(F);
    
        %vel = velocites(F);
        %Dvel = D_h*vel;
        gradR = linsolve(disc_laplacian, gradR);
       
        gradR = enforce_boudaries(F, gradR);
    
        %Compute length of gradient
        gradR_len = 0;
        for i = 1:F.num_network_index
            %WHY do I multiply by 1/2 here??
            gradR_len= gradR_len + 1/2 * norm(gradR(i, :));
        end
        fprintf("gradient length is: %d\n",gradR_len);
    
        fprintf("Energy is %d\n", J0);
        if gradR_len < 10^(-15)
            break;
        end
    
        %No backtracking loop (I guess just pure Richardson iterations).
        F = step(F, -gradR, params.alpha);
        J1 = J(F); 
    end
end


function J = J(F)
    J = 0;
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
        h = F.G(i);
    
        R1 = F.R(left, :);
        R2 = F.R(right, :);

        logR = quat_log(R1, R2);
        L_i = logR(1) * logR(2:4); %Convert to tangent vector
        
        J = J + 1/2*norm(L_i, 2)^2;
    end
end


function gradR = grad_J(F)
    gradR = zeros(F.num_network_index, 3);

    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);

        network_left  = F.network_index(left);
        network_right = F.network_index(right);

        R1 = F.R(left, :);
        R2 = F.R(right, :);

        h = F.G(i);

        logR = quat_log(R1, R2);
        L_i = logR(1) * logR(2:4); %Convert to tangent vector 

        gradR(network_left, :)  = gradR(network_left, :)  - L_i/h;
        gradR(network_right, :) = gradR(network_right, :) + L_i/h;
    end
end


function gradR = enforce_boudaries(F, gradR)
    for i = 1:size(F.fixed_nodes, 2)
        index = F.network_index(F.fixed_nodes(i));
        gradR(index, :) = zeros(1, 3);
    end
end


function F = step(F, gradR, steplength)    
    gradR = steplength*gradR;
    
    for l = 1:F.m
        index = F.network_index(l);
    
        angle = norm(gradR(index, :));
        if angle ~= 0
            axis = gradR(index, :) / angle;
            F.R(l, :) = quat_axisangle_mul(F.R(l, :), axis, angle);
        end
        %else no need to update
    end
end

function L = laplacian(F)
    L = zeros(F.num_network_index);
    
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
    
        network_left  = F.network_index(left);
        network_right = F.network_index(right);
    
        h = F.G(i);
    
        L(network_left, network_right) = L(network_left, network_right)  + 1/h^2;
        L(network_right, network_left) = L(network_right, network_left)  + 1/h^2;
        L(network_left, network_left)  = L(network_left, network_left)   - 1/h^2;
        L(network_right, network_right)= L(network_right, network_right) - 1/h^2;
    end
    
    %Without this the laplacian can be singular.
    for i = 1:size(F.fixed_nodes, 2)
        index = F.fixed_nodes(i);
        network_index = F.network_index(index);
        L(network_index, network_index) = 2 * L(network_index, network_index);
    end
    
    L = -L;
end


function vel = velocites(F)
    vel = zeros(size(F.edges, 1), 3);
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
    
        R1 = F.R(left, :);
        R2 = F.R(right, :);
    
        logR = quat_log(R1, R2);
    
        vel(i, :) = logR(1) * logR(2:4); %convert to tangent vector
    end
end


F = FiberNetwork;

F.m = 16;
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
F.R(:, 1) = 1;

% a = 0; 
% b = pi/2; 
% angle = linspace(a, b, F.m)/2;
% F.R(:, 1) = cos(angle);
% F.R(:, 4) = sin(angle);
%Hard code R
F.R(1, :) = [1, 0, 0, 0];
%F.R(2, :) = [1, 0, 0, 0];
%F.R(5, :) = [cos(pi/4), 0, sin(pi/4), 0];
F.R(F.m, :) = [cos(pi/4), 0, 0, sin(pi/4)];

draw_network(F);

J(F)

%Define solver parameters
prm = solver_param;
prm.backtracking_max = 100;
prm.beta = 0.9;
prm.alpha = 10;
prm.sigma = 1e-4;
prm.iter_max = 1000;
prm.allowed_failed_backtracks = 5;

F = laplace_solver(F, prm);

draw_network(F);