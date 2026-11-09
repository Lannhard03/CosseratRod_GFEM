%Solver for the "Kirchoff top" problem:
%For u:[0,1] \to SO(3), minimize:
%\int_[0,1] (u^T u')^T K (u^Tu') dt.
%Exact solution is just a geodesic, so using GFEM is unnecessary,
%but good test problem.

function [F, energy_hist] = kirchoff_solver(F, params)
    energy_hist = zeros(1, params.iter_max);
    %===== Outer loop =====
    J0 = J(F); 
    
    for p = 1:params.iter_max
        fprintf("Iteration %d\n", p);
        energy_hist(p) = J0;
        gradR = grad_J(F);

        %No precondition only scaling:
        %gradR = F.G(1)^2 * gradR;

        %Naive preconditioning:
        laplacian = naive_laplacian(F);
        gradR = linsolve(-laplacian, gradR);
       
        %Actual laplacian
        %laplacian = curved_laplacian(F);
        %gradR = reshape(gradR,[],1);
        %gradR = linsolve(-laplacian, gradR);
        %gradR = reshape(gradR,[], 3);

        %gradR = enforce_boudaries(F, gradR);
        
    
        %Compute length of gradient
        gradR_len = 0;
        for i = 1:F.num_network_index
            %Scale by two since we use vector norm for so(3).
            gradR_len= gradR_len + 2 * norm(gradR(i, :))^2;
        end
        
        grad_len = sqrt(gradR_len);
        fprintf("gradient length is: %d\n",grad_len);
        if grad_len < 10^(-14)
            break;
        end
    
        %If one want no backtracking loop
        F = step(F, -gradR, params.alpha);
        J0 = J(F); 

        % %===== Backtracking loop =====
        % for l = 0:params.backtracking_max
        %     steplength = params.beta^l * params.alpha;
        % 
        %     %Save previous solution
        %     R1 = F.R;
        % 
        %     %Take step in gradient direction and update energy.
        %     F = step(F, -gradR, steplength);
        %     J1 = J(F); 
        % 
        %     %Update if new energy satisfies Armijo backtracking criterion
        %     if (J0 - J1) >= params.sigma * steplength * grad_len^2
        %         break;
        %     else
        %         %Reset solution to previous one.
        %         F.R = R1;
        %     end
        % end
        % 
        % if l == params.backtracking_max
        %     fprintf("Gradient descent found no better solution iteration %d :(\n", p);
        %     failed_backtracks = failed_backtracks + 1;
        % else
        %     J0 = J1;
        %     failed_backtracks = 0;
        % end
        % if failed_backtracks > params.allowed_failed_backtracks
        %     break;
        % end

        fprintf("Energy is %d\n", J0);
    end
end


function J = J(F)
    J = 0;
    
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
    
        R1 = F.R(left, :);
        R2 = F.R(right, :);
    
        h = F.G(i);
    
        %Strain
        logR = quat_log(R1, R2);
        Rs = logR(1) * logR(2:4) / h;
    
        %Pre-strain
        Rs_pre = F.R_pre(i, :);

        %Scaled by 2 as vector I.P for so(3) is used.
        %Obs assumes K is diagonal
        curvature_term = 2*dot(F.K(i, :), (Rs - Rs_pre).^2);
    
        %Scale by 1/2 as this is done in definition of energy.
        J = J + h/2*curvature_term;
    end
end


function gradR = grad_J(F)
    gradR = zeros(F.num_network_index, 3);

    for i = 1:size(F.edges, 1)
        Rs_pre = F.R_pre(i, :);
        K = F.K(i, :);
        left = F.edges(i, 1);
        right = F.edges(i, 2);

        network_left  = F.network_index(left);
        network_right = F.network_index(right);

        R1 = F.R(left, :);
        R2 = F.R(right, :);

        h = F.G(i);

        logR = quat_log(R1, R2);
        Rs = logR(1) * logR(2:4) / h; %Convert to tangent vector
        
        %Loop over the three possible directions:
        for j = 1:3
            E = zeros(1, 3);
            E(j) = 1;

            term_1 = dot(K.* (Rs - Rs_pre),  E);

            term_2 = h * dot(K.* (Rs - Rs_pre), 1/2*cross(Rs, E));

            if network_left > 0
                gradR(network_left, j)  = gradR(network_left, j)  - term_1 + term_2;
            end
            if network_right > 0
                gradR(network_right, j) = gradR(network_right, j) + term_1 + term_2;
            end
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