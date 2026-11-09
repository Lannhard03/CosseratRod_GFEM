function F = cosserat_solver(F, params)

    J0 = trapezoid_J(F);  

    laplacian = naive_laplacian(F);
    %Assume uniform grid to test and one rod

    failed_backtracks = 0;
    %===== Outer loop =====
    for p = 1:params.iter_max
        [gradphi, gradR] = grad_J(F);
        %[gradphi, gradR] = exact_grad_J(F);
 
        %Pre-condition with discrete laplacian
        gradphi = linsolve(-laplacian, gradphi);
        %gradR = linsolve(-laplacian, gradR);

        %Pre-conditioning R with Laplacian gives worse results,
        %what is the correct way to pre-condition it?

        %curve_laplacian = curved_laplacian(F);
        %gradR = reshape(gradR, [], 1);
        %gradR = linsolve(-curve_laplacian, gradR);
        %gradR = reshape(gradR, [], 3);

    
        %Compute length of gradient
        gradphi_len = 0;
        gradR_len = 0;
        for i = 1:F.num_network_index
            gradphi_len = gradphi_len + norm(gradphi(i, :))^2;
            %GradR is in "vector form" not matrix, so we need
            %to scale the inner product by 2.
            gradR_len= gradR_len + 2 * norm(gradR(i, :))^2;
        end
        
        grad_len = sqrt(gradR_len + gradphi_len);
        if grad_len < 10^(-15)
            break;
        end
    
        %===== Backtracking loop =====
        for l = 0:params.backtracking_max
            steplength = params.beta^l * params.alpha;

            %Save previous solution
            phi1 = F.phi;
            R1 = F.R;

            %Take step in gradient direction and update energy.
            F = step(F, -gradphi, -gradR, steplength);
            J1 = trapezoid_J(F);

            %Update if new energy satisfies Armijo backtracking criterion
            if (J0 - J1) >= params.sigma * steplength * grad_len^2
                break;
            else
                %Reset solution to previous one.
                F.phi = phi1;
                F.R = R1;
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
            [gradphi, gradR] = grad_J(F);
            [gradphi_ex, gradR_ex] = exact_grad_J(F);
            %[gradphi_quad, gradR_quad] = quad_grad_J(F);
            fprintf("Relative error to exact phi grad: %d\n", norm(gradphi_ex - gradphi) / norm(gradphi_ex))
            fprintf("Relative error to exact R grad: %d\n", norm(gradR_ex - gradR) / norm(gradR))
            fprintf("Iteration %d\n", p);
            fprintf("Energy is %d\n", J0);
            J_exact = exact_J(F);
            J_quad = exact_quad_J(F);
            J_trap = trapezoid_J(F);
            %J_magnus = magnus_J(F);
            fprintf("Relative error to exact energy: %d\n", abs(J_exact - J_quad) / J_quad)
            %fprintf("Magnus difference to exact energy: %d\n", abs(J_magnus - J_quad) / J_quad)
            fprintf("Exact quad vs trapezoid energy: %d\n", abs(J_trap - J_quad) / J_quad)
            fprintf("gradient lengths are phi: %d and R: %d\n", ...
                sqrt(gradphi_len), sqrt(gradR_len));
            fprintf("Backtracking took %d iterations \n", l);
        end
    end
end


function F = step(F, gradphi, gradR, steplength)    
    gradR = steplength*gradR;

    for l = 1:F.m
        index = F.network_index(l);
        
        %Do not update fixed nodes
        if index == 0
            continue;
        end

        %Phi step
        F.phi(l, :) = F.phi(l, :) + steplength * gradphi(index, :);

        %R step
        angle = norm(gradR(index, :));
        if angle ~= 0
            axis = gradR(index, :) / angle;
            F.R(l, :) = quat_axisangle_mul(F.R(l, :), axis, angle);
        end
        %else no need to update
    end
end


function deriv = dv_Rs(j, k, h, R1, R2)
    dir = zeros(1,3);
    if k == 1
        dir(j) = 1;
    else
        dir(j) = -1;
    end
    R1_cong = R1;
    R1_cong(2:4) = -1*R1_cong(2:4);
    aux = quat2axisangle(quat_mul(R2, R1_cong));

    l = aux(2:4);
    l_p = rotate_quat(l, R1_cong);
    dir_p = rotate_quat(dir, R1);
    a = aux(1);
    %if the angle is zero the derivative is zero
    if a == 0
        cas = 1;
    else
        c = cos(a);
        s = sin(a);
        cas = c*a/s;
    end
    d = dot(l, dir_p);

    deriv = 1/(h)*(-d*l_p + cas*(d*l_p-dir) + a*cross(dir, l_p));
end

%Compute dv_Us for one element at
%the left endpoint if k = 1, and the right if k = 2.
function deriv = dv_Us(j, k, phi_p, R1, R2)
    %OBS. can change sign of dir and swap signs below to compensate!!
    dir = zeros(1,4);
    dir(j+1) = -1;

    if k == 1
        R = R1;
    elseif k == 2
        R = R2;
    end   
    R(2:4) = -R(2:4); %conjugate

    %Embed vector as imaginary quaternion
    aux1 = rotate_quat(phi_p, R);
    quat = zeros(1,4);
    quat(2:4) = aux1;

    quat = 1/2 * (quat_mul(dir,quat) - quat_mul(quat, dir));

    deriv = quat(2:4);
end

function deriv = dw_Us(j, k, h, s, R1, logR)
    dir = zeros(1,3);
    if k == 1
        dir(j) = -1/h;
    elseif k == 2
        dir(j) = 1/h;
    end

    %Take conjugate, then swap back
    R1(2:4) = -R1(2:4);
    dir = rotate_quat(dir, R1);

    %scale rotation and reverse direction
    deriv = rotate_axisangle(dir, -2*s*logR(1), logR(2:4));
end

function strain = Rs_approx(h, R1, R2)
    R1_cong = R1;
    R1_cong(2:4) = -R1_cong(2:4);
    strain_quat = quat_mul(R1_cong, R2);
    aux = quat2axisangle(strain_quat);
    strain = 1/h * aux(1) * aux(2:4);
end

function strain = Us_approx(h, phi1, phi2, R)
    phi_p = 1/h*(phi2 - phi1);
    R_cong = R;
    R_cong(2:4) = -R_cong(2:4);
    
    strain = rotate_quat(phi_p, R_cong);
end

function approx = dv_W(A, K, j, k, h, s, phi_p, R1, R2, logR, Us_pre, Rs_pre)
    %Take conjugate, then swap back
    R1(2:4) = -R1(2:4);
    dir = rotate_quat(phi_p, R1);
    R1(2:4) = -R1(2:4);
    
    %scale rotation by 2s and reverse direction
    Us = rotate_axisangle(dir, -2*s*logR(1), logR(2:4));

    factor1 = A.*(Us - Us_pre);

    %Not compute for general s in [0,1]...
    factor2 = dv_Us(j, k, phi_p, R1, R2); 
    approx = dot(factor1, factor2);
    
    %Contributions from curvature strain
    Rs_coords = Rs_approx(h, R1, R2);
    
    scaled_diff= K.*(Rs_coords - Rs_pre);
    
    dv = dv_Rs(j, k, h, R1, R2);
    
    approx = approx + 4*dot(scaled_diff, dv);
end

function approx = dw_W(A, j, k, h, s, phi_p, R1, logR, Us_pre)
    %Take conjugate, then swap back
    R1(2:4) = -R1(2:4);
    aux1 = rotate_quat(phi_p, R1);
    R1(2:4) = -R1(2:4);

    Us = rotate_axisangle(aux1, -2*s*logR(1), logR(2:4));

    factor1 = A.*(Us - Us_pre);
    
    factor2 = dw_Us(j, k, h, s, R1, logR);
    
    approx = dot(factor1, factor2);
end

function J = trapezoid_J(F)
    J = 0;
    
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);

        phi1 = F.phi(left, :);
        phi2 = F.phi(right, :);
    
        R1 = F.R(left, :);
        R2 = F.R(right, :);
    
        h = F.G(i);
    

        %Curvature term
        Rs = Rs_approx(h, R1, R2);
        Rs_pre = F.R_pre(i, :);

        %Times two from scaling of dot product
        curvature_term = 2*dot(F.K(i, :), (Rs - Rs_pre).^2);

        %1/2 half scaling from def. of energy
        J = J + h/2 * curvature_term;

        %Extension term (approximated using left and right end points)
        Us_left = Us_approx(h, phi1, phi2, R1);
        Us_right = Us_approx(h, phi1, phi2, R2);
        Us_pre = F.U_pre(i, :);
    
        extension_left  = dot(F.A(i, :), (Us_left  - Us_pre).^2);
        extension_right = dot(F.A(i, :), (Us_right - Us_pre).^2);
    
        J = J + h/2*(extension_left + extension_right)/2;

    end
end

function J = magnus_J(F)
    J = 0;
    
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
    
        phi1 = F.phi(left, :);
        phi2 = F.phi(right, :);
    
        R1 = F.R(left, :);
        R2 = F.R(right, :);
    
        h = F.G(i);
    
    
        %Curvature term
        Rs = Rs_approx(h, R1, R2);
        Rs_pre = F.R_pre(i, :);
    
        %Times two from scaling of dot product
        curvature_term = 2*dot(F.K(i, :), (Rs - Rs_pre).^2);
    
        %1/2 half scaling from def. of energy
        J = J + h/2 * curvature_term;
    
        %Extension term (approximated using left and right end points)
        Us = Us_approx(h, phi1, phi2, R1);
        Us_pre = F.U_pre(i, :);
        diff = Us-Us_pre;

        %Matrix times diff done fast using diagonal structure
        %Magnus approximation of integral.
        w = F.A(i, :) .* diff;
        y = cross(Rs, diff);
        x = dot(diff, w) + h^2/3 * dot(F.A(i, :), y.^2); 
    
        J = J + h/2* x;
    end
end

function J = exact_quad_J(F)
    J = 0;
    
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
    
        phi1 = F.phi(left, :);
        phi2 = F.phi(right, :);
    
        R1 = F.R(left, :);
        R2 = F.R(right, :);
    
        h = F.G(i);
    
        Rs = Rs_approx(h, R1, R2);
    
        Us_pre = F.U_pre(i, :);
        Rs_pre = F.R_pre(i, :);
    

        %Extension term (integral computed exactly...)
        phi_p = (phi2 - phi1) / h;

        logR = quat_log(R1, R2);
        
        n_i = logR(2:4);
        theta_i = logR(1);
        
        R1(2:4) = -R1(2:4);
        Us_i = rotate_quat(phi_p, R1); %Strain at left endpoint
        R1(2:4) = -R1(2:4);

        %diff = Us_i - Us_pre;

        %Quadrature! 
        integrand = @(x) F.A(i, :) * (rotate_axisangle(Us_i, -theta_i*x, n_i) - Us_pre).^2.';
        
        extension_term = integral(integrand, 0, 1); %, RelTol=1e-2
        
        %Curvature term
        %*2 from scaling of I.P, and then scale 1/2 from def. of energy.
        curvature_term = 2*dot(F.K(i, :), (Rs - Rs_pre).^2);
        J = J + h/2*curvature_term + h/2*extension_term;
    end
end

function y = mysinc(x)
    y = ones(size(x));
    nz = (x ~= 0);
    y(nz) = sin(x(nz)) ./ (x(nz));
end

function J = exact_J(F)
    J = 0;
    
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
    
        phi1 = F.phi(left, :);
        phi2 = F.phi(right, :);
    
        R1 = F.R(left, :);
        R2 = F.R(right, :);
    
        h = F.G(i);
    
        Rs = Rs_approx(h, R1, R2);
    
        Us_pre = F.U_pre(i, :);
        Rs_pre = F.R_pre(i, :);
    
    
        %Extension term
        phi_p = (phi2 - phi1) / h;
    
        logR = quat_log(R1, R2);
    
        n = logR(2:4);
        theta = logR(1);
   
        R1(2:4) = -R1(2:4);
        Us_i = rotate_quat(phi_p, R1); %Strain at left endpoint
        R1(2:4) = -R1(2:4);
    
        %Us_i = Us_i - Us_pre;
        tang = dot(Us_i, n)*n;
        orth1 = Us_i - tang;
        orth2 = cross(n, orth1);
        
        extension_term = 0;
        %Pre-strain term
        w1 = F.A(i, :) .* Us_pre;
        extension_term = extension_term + dot(w1, Us_pre);

        %Mixed terms
        extension_term = extension_term - 2 * dot(w1, tang);
        extension_term = extension_term - 2 * mysinc(theta) * dot(w1, orth1);
        extension_term = extension_term + 2 * theta/2 * mysinc(theta/2)^2 * dot(w1, orth2);

        %Strain terms
        extension_term = extension_term + dot(F.A(i, :) .* tang, tang);
        extension_term = extension_term + 2 * mysinc(theta) * dot(F.A(i, :) .* tang, orth1); 
        extension_term = extension_term - theta*mysinc(theta/2)^2 * dot(F.A(i, :) .* tang, orth2); 
        extension_term = extension_term + (1/2 + mysinc(2*theta)/2) * dot(F.A(i, :) .* orth1, orth1);
        extension_term = extension_term - theta*mysinc(theta)^2 * dot(F.A(i, :) .* orth2, orth1);
        extension_term = extension_term + (1/2 - mysinc(2*theta)/2) * dot(F.A(i, :) .* orth2, orth2);

        %Curvature term
        %*2 from scaling of I.P, and then scale 1/2 from def. of energy.
        curvature_term = 2*dot(F.K(i, :), (Rs - Rs_pre).^2);
        J = J + h/2*curvature_term + h/2*extension_term;
    end
end

function [gradphi, gradR] = grad_J(F)
    gradphi = zeros(F.num_network_index, 3);
    gradR = zeros(F.num_network_index, 3);
    
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);

        network_left  = F.network_index(left);
        network_right = F.network_index(right);

        phi1 = F.phi(left, :);
        phi2 = F.phi(right, :);

        R1 = F.R(left, :);
        R2 = F.R(right, :);

        h = F.G(i);

        Us_pre = F.U_pre(i, :);
        Rs_pre = F.R_pre(i, :);

        phi_p = (phi2 - phi1) / h;
    
        logR = quat_log(R1, R2);

        for j = 1:3
            %===dw_W derivatives
            %Compute change at the left and right endpoints, 
            %when w changes in j direction at left endpoint
            if network_left > 0
                dw_left_left  = dw_W(F.A(i, :), j, 1, h, 0, phi_p, R1, logR, Us_pre);
                dw_right_left = dw_W(F.A(i, :), j, 1, h, 1, phi_p, R1, logR, Us_pre);

                dw_total = h/2*(dw_left_left + dw_right_left);
                gradphi(network_left, j) = gradphi(network_left, j) + dw_total; 
            end

                
            %Compute change at the left and right endpoints, 
            %when w changes in j direction at right endpoint
            if network_right > 0
                dw_left_right  = dw_W(F.A(i, :), j, 2, h, 0, phi_p, R1, logR, Us_pre);
                dw_right_right = dw_W(F.A(i, :), j, 2, h, 1, phi_p, R1, logR, Us_pre);

                dw_total = h/2*(dw_left_right + dw_right_right);
                gradphi(network_right, j) = gradphi(network_right, j) + dw_total; 
            end


            %===dv_W derivatives
            if network_left > 0
                dv_left_left  = dv_W(F.A(i, :), F.K(i, :), j, 1, h, 0, phi_p, R1, R2, logR, Us_pre, Rs_pre);
                dv_right_left = dv_W(F.A(i, :), F.K(i, :), j, 1, h, 1, phi_p, R1, R2, logR, Us_pre, Rs_pre);
    
                total_dv = h/2*(dv_left_left + dv_right_left);
                gradR(network_left, j) = gradR(network_left, j) + total_dv;
            end

            if network_right > 0
                dv_left_right  = dv_W(F.A(i, :), F.K(i, :), j, 2, h, 0, phi_p, R1, R2, logR, Us_pre, Rs_pre);
                dv_right_right = dv_W(F.A(i, :), F.K(i, :), j, 2, h, 1, phi_p, R1, R2, logR, Us_pre, Rs_pre);
    
                total_dv = h/2*(dv_left_right + dv_right_right);
                gradR(network_right, j) = gradR(network_right, j) + total_dv;
            end
        end
    end 
end

function [gradphi, gradR] = quad_grad_J(F)
    gradphi = zeros(F.num_network_index, 3);
    gradR = zeros(F.num_network_index, 3);
    
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
    
        K = F.K(i, :);
        A = F.A(i, :);
    
        network_left  = F.network_index(left);
        network_right = F.network_index(right);
    
        phi1 = F.phi(left, :);
        phi2 = F.phi(right, :);
    
        R1 = F.R(left, :);
        R2 = F.R(right, :);
    
        h = F.G(i);
    
        Us_pre = F.U_pre(i, :);
        Rs_pre = F.R_pre(i, :);
    
        phi_p = (phi2 - phi1) / h;
    
        logR = quat_log(R1, R2);
    
        Rs = logR(1) * logR(2:4) / h;
    
        n = logR(2:4);
        theta = logR(1);
        v = theta*n;
    
        R1(2:4) = -R1(2:4);
        Us_i = rotate_quat(phi_p, R1); %Strain at left endpoint
        R1(2:4) = -R1(2:4);

        for j = 1:3
            E = zeros(1, 3);
            E(j) = 1;
    
            %=========gradR==============
    
            %Contribution from curvature strain:
            term_1 = dot(K.* (Rs - Rs_pre),  E);
            term_2 = h * dot(K.* (Rs - Rs_pre), 1/2*cross(Rs, E));
    
            %Contribution from extension strain:
            %Quadrature! 

            hat = @(w) [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
            eps = 1e-6;
            M = diag(A);
            

            if network_left > 0
                delta = -E + 0.5*cross(v, E);
                f = @(t) ((expm(-t*hat(v+eps*delta))*Us_i.' - expm(-t*hat(v-eps*delta))*Us_i.')/(2*eps))' ...
                    * M * (expm(-t*hat(v))*Us_i.' - Us_pre.');

                I_num = integral(@(t) arrayfun(f,t), 0, 1, 'AbsTol',1e-2,'RelTol',1e-2);

                gradR(network_left, j)  = gradR(network_left, j)  - term_1 + term_2 + I_num;
            end
            if network_right > 0
                delta = E + 0.5*cross(v, E);
                f = @(t) ((expm(-t*hat(v+eps*delta))*Us_i.' - expm(-t*hat(v-eps*delta))*Us_i.')/(2*eps))' ...
                    * M * (expm(-t*hat(v))*Us_i.' - Us_pre.');

                I_num = integral(@(t) arrayfun(f,t), 0, 1, 'AbsTol',1e-2,'RelTol',1e-2);

                gradR(network_right, j) = gradR(network_right, j) + term_1 + term_2 + I_num;
            end
    
            %=========gradphi==============
    
            R1(2:4) = -R1(2:4);
            gamma = rotate_quat(E, R1); %Strain at left endpoint
            R1(2:4) = -R1(2:4);
    
            %Quadrature! 
            integrand = @(x) dot((rotate_axisangle(gamma, -theta*x, n)), (F.A(i, :).*(rotate_axisangle(Us_i, -theta*x, n) - Us_pre)));

            extension_term = integral(integrand, 0, 1, Vectorized=0); %, RelTol=1e-2
    
            %Contribution from extension strain:
            if network_left > 0
                gradphi(network_left, j)  = gradphi(network_left, j)  - extension_term;
            end
            if network_right > 0
                gradphi(network_right, j) = gradphi(network_right, j) + extension_term;
            end
        end
    end 
end

function [gradphi, gradR] = exact_grad_J(F)
    gradphi = zeros(F.num_network_index, 3);
    gradR = zeros(F.num_network_index, 3);
    
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);

        K = F.K(i, :);
        A = F.A(i, :);

        network_left  = F.network_index(left);
        network_right = F.network_index(right);

        phi1 = F.phi(left, :);
        phi2 = F.phi(right, :);

        R1 = F.R(left, :);
        R2 = F.R(right, :);

        h = F.G(i);

        Us_pre = F.U_pre(i, :);
        Rs_pre = F.R_pre(i, :);

        phi_p = (phi2 - phi1) / h;

        logR = quat_log(R1, R2);

        Rs = logR(1) * logR(2:4) / h;
        v = logR(1) * logR(2:4);

        n = logR(2:4);
        theta = logR(1);

        R1(2:4) = -R1(2:4);
        Us_i = rotate_quat(phi_p, R1); %Strain at left endpoint
        R1(2:4) = -R1(2:4);

        for j = 1:3
            E = zeros(1, 3);
            Z = zeros(1, 3);
            E(j) = 1;

            %=========gradR==============

            %Contribution from curvature strain:
            term_1 = 2*dot(K.* (Rs - Rs_pre),  E);
            term_2 = 2*h * dot(K.* (Rs - Rs_pre), 1/2*cross(Rs, E));
            if network_left > 0
                gradR(network_left, j)  = gradR(network_left, j)  - term_1 + term_2;
            end
            if network_right > 0
                gradR(network_right, j) = gradR(network_right, j) + term_1 + term_2;
            end

            if network_left > 0
                left = rot_integral_variation(v, Z, E, Us_i, Us_pre, diag(A));

                gradR(network_left, j)  = gradR(network_left, j)  + 2*h*left;
            end
            if network_right > 0
                right = rot_integral_variation(v, E, Z, Us_i, Us_pre, diag(A));

                gradR(network_right, j) = gradR(network_right, j) + 2*h*right;
            end
            
            %=========gradphi==============
            
            R1(2:4) = -R1(2:4);
            gamma = rotate_quat(E, R1); %Strain at left endpoint
            R1(2:4) = -R1(2:4);

            extension_term = 0;

            %Mixed terms
            extension_term = extension_term - I2(-theta, n, gamma, Us_pre, A);

            %Strain terms
            extension_term = extension_term + I1(-theta, n, gamma, Us_i, A);

            %Contribution from extension strain:
            if network_left > 0
                gradphi(network_left, j)  = gradphi(network_left, j)  - extension_term;
            end
            if network_right > 0
                gradphi(network_right, j) = gradphi(network_right, j) + extension_term;
            end
        end
    end 
end

function result = I2(theta, n, x, y, A)
    tang = dot(x, n)*n;
    orth1 = x - tang;
    orth2 = cross(n, orth1);
    
    result = 0;
    w1 = A .* y;
    
    %Mixed terms
    result = result + dot(w1, tang);
    result = result + mysinc(theta) * dot(w1, orth1);
    result = result + theta/2 * mysinc(theta/2)^2 * dot(w1, orth2);
end

function result = I1(theta, n, x, y, A)
    tangx = dot(x, n)*n;
    orth1x = x - tangx;
    orth2x = cross(n, orth1x);

    tangy = dot(y, n)*n;
    orth1y = y - tangy;
    orth2y = cross(n, orth1y);
    
    result = 0;
    result = result + dot(A.*tangx, tangy);
    result = result + mysinc(theta) * (dot(A.*tangx, orth1y) + dot(A.*tangy, orth1x));
    result = result + theta/2 * mysinc(theta/2)^2 * (dot(A.*tangx, orth2y) + dot(A.*tangy, orth2x));
    result = result + (1/2 + mysinc(2*theta)/2) * dot(A.*orth1x, orth1y);
    result = result + theta/2*mysinc(theta)^2 * (dot(A.*orth1x, orth2y) + dot(A.*orth1y, orth2x));
    result = result + (1/2 - mysinc(2*theta)/2) * dot(A.*orth2x, orth2y);
end

%Written by Claude: (SLOW AS IT CONVERTS TO MATRIX)
function [x1,x2,x3,x4,x5] = rot_deriv_vectors(v, w_next, w_cur, xi)
    % Vectors such that
    %  d/ds|_{s=0} [ expm(-t*hat(x(s))) * xi(s) ],  xi(s) = expm(-s*hat(w_cur))*xi,
    %    = (1-cos(th t))x1 - t sin(th t)x2 + t cos(th t)x3 + sin(th t)x4 + cos(th t)x5
    % with x(0)=v, x'(0) = w_next - w_cur + 0.5*cross(v, w_next + w_cur).
    v=v(:); xi=xi(:); w_next=w_next(:); w_cur=w_cur(:);
    th    = norm(v);
    delta = w_next - w_cur + 0.5*cross(v, w_next + w_cur);
    eta   = cross(xi, w_cur);               % d/ds xi(s) at s=0
    l     = v/th;
    ld    = l'*delta;   lx = l'*xi;
    dperp = delta - ld*l;
    xperp = xi - lx*l;
    
    % axis-variation part
    x1 = ( lx*dperp + (dperp'*xi)*l )/th;
    x2 = ld*xperp;
    x3 = -ld*cross(l, xi);
    x4 = ( lx*cross(l, delta) - (l'*cross(delta, xi))*l )/th;
    
    % xi(s) part: R*eta = (1-cos)(l.eta) l + cos*eta - sin*(l x eta)
    x1 = x1 + (l'*eta)*l;
    x4 = x4 - cross(l, eta);
    x5 = eta;
end


function I = rot_integral_four(x1, x2, x3, x4, x5, M, v, xi, xit)
    % Exact value of
    %  int_0^1 ((1-cos(th t))x1 - t sin(th t)x2 + t cos(th t)x3 + sin(th t)x4
    %           + cos(th t)x5)' * M * (expm(-t*hat(v))*xi - xit) dt,   th = norm(v).
    x1=x1(:); x2=x2(:); x3=x3(:); x4=x4(:); x5=x5(:);
    xi=xi(:); xit=xit(:); v=v(:);
    th = norm(v);
    if th < 1e-12, l = [0;0;1]; else, l = v/th; end
    K = [0 -l(3) l(2); l(3) 0 -l(1); -l(2) l(1) 0];
    P = l*l';  Q = eye(3) - P;
    
    if th < 1e-2
        c1 = 1 - th^2/6 + th^4/120;
        s1 = th/2 - th^3/24 + th^5/720;
        c2 = 1 - th^2/3 + th^4/15;
        s2 = th^2/3 - th^4/15;
        sc = th/2 - th^3/6 + th^5/45;
        Ts = th/3 - th^3/30 + th^5/840;
        Tc = 1/2 - th^2/8 + th^4/144;
        Tsc = th/3 - 2*th^3/15 + 2*th^5/105;
        Tcc = 1/2 - th^2/4 + th^4/18;
        Tss = th^2/4 - th^4/18;
    else
        c1 = sin(th)/th;
        s1 = (1-cos(th))/th;
        c2 = 1/2 + sin(2*th)/(4*th);
        s2 = 1/2 - sin(2*th)/(4*th);
        sc = sin(th)^2/(2*th);
        Ts = (sin(th) - th*cos(th))/th^2;
        Tc = (cos(th) + th*sin(th) - 1)/th^2;
        Tsc = (sin(2*th) - 2*th*cos(2*th))/(8*th^2);
        E = (cos(2*th) + 2*th*sin(2*th) - 1)/(8*th^2);
        Tcc = 1/4 + E;
        Tss = 1/4 - E;
    end
    
    X = [x1 x2 x3 x4 x5];
    A = X'*M*(P*xi - xit);
    q = X'*M*(Q*xi);
    r = X'*M*(K*xi);
    
    I =   A(1)*(1-c1) + q(1)*(c1-c2) - r(1)*(s1-sc) ...
        - A(2)*Ts     - q(2)*Tsc     + r(2)*Tss ...
        + A(3)*Tc     + q(3)*Tcc     - r(3)*Tsc ...
        + A(4)*s1     + q(4)*sc      - r(4)*s2 ...
        + A(5)*c1     + q(5)*c2      - r(5)*sc;
end

function dJ = rot_integral_variation(v, w_next, w_cur, xi, xit, M)
    % d/ds|_{s=0} of  1/2 int_0^1 (E(s,t) xi(s) - xit)' M (E(s,t) xi(s) - xit) dt
    % (M symmetric; the 1/2 cancels the factor 2 from the product rule.)
    v=v(:); xi=xi(:); xit=xit(:); w_next=w_next(:); w_cur=w_cur(:);
    th = norm(v);
    if th < 1e-6
        delta = w_next - w_cur + 0.5*cross(v, w_next + w_cur);
        eta   = cross(xi, w_cur);
        dJ = (eta - 0.5*cross(delta, xi))' * M * (xi - xit);
        return
    end
    [x1,x2,x3,x4,x5] = rot_deriv_vectors(v, w_next, w_cur, xi);
    dJ = rot_integral_four(x1,x2,x3,x4,x5, M, v, xi, xit);
end