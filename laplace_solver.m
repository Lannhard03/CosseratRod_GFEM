%Solver for the "dirichlet energy" problem:
%For u:[0,1] \to SO(3), minimize:
%\int_[0,1] <u', u'> dt.
%Exact solution is just a geodesic, so using GFEM is unnecessary,
%but good test problem.

function [F, energy_hist] = laplace_solver(F, params)
    energy_hist = zeros(1, params.iter_max);
    %===== Outer loop =====
    J1 = J(F); 
    
    for p = 1:params.iter_max
        fprintf("Iteration %d\n", p);
        energy_hist(p) = J1;
        gradR = grad_J(F);

        %No precondition only scaling:
        %gradR = F.G(1)^2 * gradR;

        %Naive preconditioning:
        %laplacian = naive_laplacian(F);
        %gradR = linsolve(-laplacian, gradR);
       
        %Actual laplacian
        laplacian = curved_laplacian(F);
        gradR = reshape(gradR,[],1);
        gradR = linsolve(-laplacian, gradR);
        gradR = reshape(gradR,[], 3);



        %gradR = enforce_boudaries(F, gradR);
        
    
        %Compute length of gradient
        gradR_len = 0;
        for i = 1:F.num_network_index
            %WHY do I multiply by 1/2 here??
            gradR_len= gradR_len + 1/2 * norm(gradR(i, :));
        end
        fprintf("gradient length is: %d\n",gradR_len);
    
        
        if gradR_len < 10^(-14)
            break;
        end
    
        %No backtracking loop (I guess just pure Richardson iterations).
        F = step(F, -gradR, params.alpha);
        J1 = J(F); 

        fprintf("Energy is %d\n", J1);
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
        L_i = logR(1) * logR(2:4) / h; %Convert to tangent vector
        
        J = J + h*dot(L_i, L_i);
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
        
        if network_left > 0
            gradR(network_left, :)  = gradR(network_left, :)  - L_i/h;
        end
        if network_right > 0
            gradR(network_right, :) = gradR(network_right, :) + L_i/h;
        end
    end
end


function F = step(F, gradR, steplength)    
    gradR = steplength*gradR;
    
    for l = 1:F.m
        index = F.network_index(l);
        %Do not update fixed nodes
        if index == 0
            continue;
        end
    
        angle = norm(gradR(index, :));
        if angle ~= 0
            axis = gradR(index, :) / angle;
            F.R(l, :) = quat_axisangle_mul(F.R(l, :), axis, angle);
        end
        %else no need to update
    end
end

