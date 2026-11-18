function F = thimoshenko_solver(F, params)
    J0 = J(F);  

    laplacian = naive_laplacian(F);
    
    failed_backtracks = 0;
    %===== Outer loop =====
    for p = 1:params.iter_max

        %F = F.update_network();

        [gradu, gradr] = grad_J(F);
 
        %Pre-condition with discrete laplacian
        gradu = linsolve(-laplacian, gradu);
        gradr = linsolve(-laplacian, gradr);
    
        %Compute length of gradient
        gradu_len = 0;
        gradr_len = 0;
        for i = 1:F.num_network_index
            gradu_len = gradu_len + norm(gradu(i, :))^2;
            %GradR is in "vector form" not matrix, so we need
            %to scale the inner product by 2.
            gradr_len= gradr_len + norm(gradr(i, :))^2;
        end
        
        grad_len = sqrt(gradr_len + gradu_len);
        if grad_len < 10^(-15)
            break;
        end
    
        %===== Backtracking loop =====
        for l = 0:params.backtracking_max
            steplength = params.beta^l * params.alpha;

            %Save previous solution
            u1 = F.u;
            r1 = F.r;

            %Take step in gradient direction and update energy.
            F = step(F, -gradu, -gradr, steplength);
            J1 = J(F);

            %Update when new energy satisfies Armijo backtracking criterion
            if (J0 - J1) >= params.sigma * steplength * grad_len^2
                break;
            else
                %Reset solution to previous one.
                F.u = u1;
                F.r = r1;
            end
        end
    
        if l == params.backtracking_max
            fprintf("Gradient descent found no better solution iteration %d :(\n", p);
            failed_backtracks = failed_backtracks + 1;
        else
            J0 = J1;
            
            failed_backtracks = 0;
        end
        if failed_backtracks > params.allowed_failed_backtracks
            break;
        end
        
        %Print info every 10 iterations
        if mod(p, 10) == 0
            fprintf("Iteration %d\n", p);
            fprintf("Energy is %d\n", J0);
            fprintf("gradient lengths are phi: %d and R: %d\n", ...
                sqrt(gradu_len), sqrt(gradr_len));
            fprintf("Backtracking took %d iterations \n", l);
        end
    end
end


function F = step(F, gradu, gradr, steplength)    
    for l = 1:F.m
        index = F.network_index(l);
        
        %Do not update fixed nodes
        if index == 0
            continue;
        end

        %u step
        F.u(l, :) = F.u(l, :) + steplength * gradu(index, :);

        %r step
        F.r(l, :) = F.r(l, :) + steplength * gradr(index, :);
    end
end

function J = J(F)
    left = F.edges(:, 1);
    right = F.edges(:, 2);

    r_l = F.r(left, :);
    r_r = F.r(right, :);

    u_l = F.u(left, :);
    u_r = F.u(right, :);

    %Derivatives: u prime and r prime
    u_p = (u_r - u_l) ./ F.G;
    r_p = (r_r - r_l) ./ F.G;

    r_cross_l = cross(F.n , r_l, 2);
    r_cross_r = cross(F.n , r_r, 2);


    %Compute each term
    w1 = F.A .* u_p;

    T1 = sum(F.G .* dot(w1, u_p, 2));
    T2 = sum(F.G .* dot(w1, r_cross_r + r_cross_l, 2));

    T3 = sum(F.G/3 .* dot(F.A .* r_cross_l, r_cross_l, 2));
    T4 = sum(F.G/3 .* dot(F.A .* r_cross_r, r_cross_r, 2));
    T5 = sum(F.G/3 .* dot(F.A .* r_cross_l, r_cross_r, 2));

    T6 = sum(F.G .* dot(F.K .* r_p, r_p, 2));

    J = (T1 + T2 + T3 + T4 + T5 + T6) / 2;
end

function [gradu, gradr] = grad_J(F)
    gradu = zeros(F.num_network_index, 3);
    gradr = zeros(F.num_network_index, 3);
    
    left = F.edges(:, 1);
    right = F.edges(:, 2);

    r_l = F.r(left, :);
    r_r = F.r(right, :);

    u_l = F.u(left, :);
    u_r = F.u(right, :);

    %Derivatives: u prime and r prime
    u_p = (u_r - u_l) ./ F.G;
    u_p_cross = cross(F.A .* u_p, F.n, 2);

    r_p = (r_r - r_l) ./ F.G;

    r_cross_l = cross(F.n , r_l, 2);
    r_cross_r = cross(F.n , r_r, 2);
    r_cross_ave = (r_cross_r + r_cross_l)/2;

    r_2cross_l = cross(F.A .* r_cross_l, F.n, 2);
    r_2cross_r = cross(F.A .* r_cross_r, F.n, 2);

    %Terms for u gradient
    gradu_terms = F.A .* (u_p + r_cross_ave);

    %Terms for r gradient, split into terms which depend on the nodal
    %value at the left or right endpoint.
    gradr_left  = F.G.* (u_p_cross/2 + r_2cross_l / 3 + r_2cross_r / 6 - F.K .* r_p);

    gradr_right = F.G.* (u_p_cross/2 + r_2cross_r / 3 + r_2cross_l / 6 + F.K .* r_p);

    for i = 1:size(F.edges, 1)
        network_left  = F.network_index(left(i));
        network_right = F.network_index(right(i));

        if network_left > 0
            gradu(network_left, :) = gradu(network_left, :) - gradu_terms(i, :);
            gradr(network_left, :) = gradr(network_left, :) + gradr_left(i, :);
        end
        if network_right > 0
            gradu(network_right, :) = gradu(network_right, :) + gradu_terms(i, :);
            gradr(network_right, :) = gradr(network_right, :) + gradr_right(i, :);
        end
    end
end
