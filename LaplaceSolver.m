%Solver for the "harmonic energy" problem:
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
        naive_disc_laplacian = naive_laplacian(F);
        gradR = linsolve(naive_disc_laplacian, gradR);
       
        %Actual laplacian
        %disc_laplacian = laplacian(F);
        %gradR = reshape(gradR,[],1);
        %gradR = linsolve(disc_laplacian, gradR);
        %gradR = reshape(gradR,[], 3);



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
        
        J = J + h*norm(L_i, 2)^2;
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

function L = laplacian(F)
    %Need an element for each direction.
    n_net = F.num_network_index;
    L = zeros(3 * n_net);

    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
    
        n_left  = F.network_index(left);
        n_right = F.network_index(right);
    
        h = F.G(i);

        R1 = F.R(left, :);
        R2 = F.R(right, :);
        
        logR = quat_log(R1, R2);
        L_i = logR(1) * logR(2:4); %Convert to tangent vector 
        
        L_norm = norm(L_i);
        if L_norm == 0
            %Limit case of the below recovers naive laplacian.
            for j = 0:2
                if n_left == 0 && n_right == 0
                    continue;
                elseif n_left == 0 %In this case network_right != 0
                     L(n_right + n_net*j, n_right + n_net*j) = L(n_right + n_net*j, n_right + n_net*j)  - 1/h;
                elseif n_right == 0 %in this case network_left != 0
                    L(n_left  + n_net*j, n_left  + n_net*j) = L(n_left  + n_net*j, n_left  + n_net*j)  - 1/h;
                else
                    L(n_right + n_net*j, n_right + n_net*j) = L(n_right + n_net*j, n_right + n_net*j)  - 1/h;
                    L(n_left  + n_net*j, n_left  + n_net*j) = L(n_left  + n_net*j, n_left  + n_net*j)  - 1/h;

                    L(n_left  + n_net*j, n_right + n_net*j) = L(n_left  + n_net*j, n_right + n_net*j)  + 1/h;
                    L(n_right + n_net*j, n_left  + n_net*j) = L(n_right + n_net*j, n_left  + n_net*j)  + 1/h;
                end
            end
        else 

            K = eye(3) - ((L_i/h).'*(L_i/h))/(L_norm/h)^2;
            S = (L_i/h).'*(L_i/h)/(L_norm/h)^2;
    
            b = L_norm/(2*sqrt(2));
            s = sin(b);
            c = cos(b);
    
            main =  K/(h*s^2) * b/2 * (b + s*c) + S/h;
            side = -K/(h*s^2) * b/2 * (s + b*c) - S/h;
            
            
            for j = 0:2
                for l = 0:2
                    if n_left == 0 && n_right == 0
                        continue;
                    elseif n_left == 0 %In this case network_right != 0
                        L(n_right + n_net*j, n_right + n_net*j) = L(n_right + n_net*j, n_right + n_net*j)  + main(j+1, l+1);
                    elseif n_right == 0 %in this case network_left != 0
                        L(n_left  + n_net*j, n_left  + n_net*j) = L(n_left  + n_net*j, n_left  + n_net*j)  + main(j+1, l+1);
                    else
                        L(n_right + n_net*j, n_right + n_net*j) = L(n_right + n_net*j, n_right + n_net*j)  + main(j+1, l+1);
                        L(n_left  + n_net*j, n_left  + n_net*j) = L(n_left  + n_net*j, n_left  + n_net*j)  + main(j+1, l+1);

                        L(n_left  + n_net*j, n_right + n_net*j) = L(n_left  + n_net*j, n_right + n_net*j)  + side(j+1, l+1);
                        L(n_right + n_net*j, n_left  + n_net*j) = L(n_right + n_net*j, n_left  + n_net*j)  + side(j+1, l+1);
                    end
                end
            end
        end
    end
end


function L = naive_laplacian(F)
    L = zeros(F.num_network_index);
    
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
    
        network_left  = F.network_index(left);
        network_right = F.network_index(right);
    
        h = F.G(i);
    
        if network_left == 0 && network_right == 0
            continue;
        elseif network_left == 0 %In this case network_right != 0
            L(network_right, network_right)= L(network_right, network_right) - 1/h;
        elseif network_right == 0 %in this case network_left != 0
            L(network_left, network_left)  = L(network_left, network_left)   - 1/h;
        else
            L(network_left, network_left)  = L(network_left, network_left)   - 1/h;
            L(network_right, network_right)= L(network_right, network_right) - 1/h;

            L(network_left, network_right) = L(network_left, network_right)  + 1/h;
            L(network_right, network_left) = L(network_right, network_left)  + 1/h;
        end
    end
    
    L = -L;
end