function F = network_solver(F, params)

    J0 = J(F);    

    discrete_laplacian = disc_laplacian(F);

    failed_backtracks = 0;
    %===== Outer loop =====
    for p = 1:params.iter_max
        fprintf("Iteration %d\n", p);

        [gradphi, gradR] = grad_J(F);
 
        %Pre-condition with discrete laplacian
        gradphi = linsolve(discrete_laplacian, gradphi);

        %Pre-conditioning R with Laplacian gives worse results,
        %what is the correct way to pre-condition it?
        %discrete_quat_laplacian = disc_quat_laplacian(h, m_tot, rods(1).R);
        %gradR_col = reshape(gradR.', [], 1);
        %gradR_col = linsolve(discrete_quat_laplacian, gradR_col);
        %gradR = reshape(gradR_col, 3, []).';

        [gradphi,gradR] = enforce_boudaries(F,gradphi,gradR);
    
        %Compute length of gradient
        gradphi_len = 0;
        gradR_len = 0;
        for i = 1:F.num_network_index
            gradphi_len = gradphi_len + norm(gradphi(i, :));
            %WHY do I multiply by 1/2 here??
            gradR_len= gradR_len + 1/2 * norm(gradR(i, :));
        end
        fprintf("gradient lengths are phi: %d and R: %d\n", gradphi_len, gradR_len);
  
        fprintf("Energy is %d\n", J0);
        if gradR_len + gradphi_len < 10^(-15)
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
            J1 = J(F); 

            %Update if new energy is lower
            if (J0 - J1) >= params.sigma*steplength*(gradphi_len + gradR_len)
                break;
            else
                %Reset solution to previous one.
                F.phi = phi1;
                F.R = R1;
            end
        end
    
        fprintf("Backtracking took %d iterations \n", l);
    
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
    end
end

function F = step(F, gradphi, gradR, steplength)    
    gradR = steplength*gradR;

    for l = 1:F.m
        index = F.network_index(l);
        %Phi step
        F.phi(l, :) = F.phi(l, :) + steplength * gradphi(index, :);

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
    R1(2:4) = -R1(2:4);

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
    
    approx = approx + dot(scaled_diff, dv);
end

function approx = dw_W(A, j, k, h, s, phi_p, R1, logR, Us_pre)
    %Take conjugate, then swap back
    R1(2:4) = -R1(2:4);
    aux = rotate_quat(phi_p, R1);
    R1(2:4) = -R1(2:4);

    Us = rotate_axisangle(aux, -2*s*logR(1), logR(2:4));

    factor1 = A.*(Us - Us_pre);
    
    factor2 = dw_Us(j, k, h, s, R1, logR);
    
    approx = dot(factor1, factor2);
end

function J = J(F)
    J = 0;
    
    %Use Trapezoid rule to sum energy contributions from each edge in
    %network
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
        phi1 = F.phi(left, :);
        phi2 = F.phi(right, :);
    
        R1 = F.R(left, :);
        R2 = F.R(right, :);
    
        h = F.G(i);
    
        Rs = Rs_approx(h, R1, R2);
    
        %Left endpoint
        Us_pre = F.Us_pre(i, :);
        Rs_pre = F.Rs_pre(i, :);
    
        Us = Us_approx(h, phi1, phi2, R1);
    
        extension_term = dot(F.A(i, :), (Us - Us_pre).^2);
        %Why times 4 here??
        curvature_term = 4*dot(F.K(i, :), (Rs - Rs_pre).^2);
    
        J = J + h/2*(extension_term + curvature_term);
    
        %Right endpoints
        Us = Us_approx(h, phi1, phi2, R2);
    
        extension_term = dot(F.A(i, :), (Us - Us_pre).^2);
        %Why times 4 here??
        curvature_term = 4*dot(F.K(i, :), (Rs - Rs_pre).^2);

        J = J + h/2*(extension_term + curvature_term);
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

        Us_pre = F.Us_pre(i, :);
        Rs_pre = F.Rs_pre(i, :);

        phi_p = (phi2 - phi1) / h;
    
        logR = quat_log(R1, R2);

        for j = 1:3
            %===dw_W derivatives
            %Compute change at the left and right endpoints, 
            %when w changes in j direction at left endpoint
            dw_left_left  = dw_W(F.A(i, :), j, 1, h, 0, phi_p, R1, logR, Us_pre);
            dw_right_left = dw_W(F.A(i, :), j, 1, h, 1, phi_p, R1, logR, Us_pre);

            dw_total = h/2*(dw_left_left + dw_right_left);
            gradphi(network_left, j) = gradphi(network_left, j) + dw_total; 
                
            %Compute change at the left and right endpoints, 
            %when w changes in j direction at right endpoint
            dw_left_right  = dw_W(F.A(i, :), j, 2, h, 0, phi_p, R1, logR, Us_pre);
            dw_right_right = dw_W(F.A(i, :), j, 2, h, 1, phi_p, R1, logR, Us_pre);

            dw_total = h/2*(dw_left_right + dw_right_right);
            gradphi(network_right, j) = gradphi(network_right, j) + dw_total; 

            %===dv_W derivatives
            dv_left_left  = dv_W(F.A(i, :), F.K(i, :), j, 1, h, 0, phi_p, R1, R2, logR, Us_pre, Rs_pre);
            dv_right_left = dv_W(F.A(i, :), F.K(i, :), j, 1, h, 1, phi_p, R1, R2, logR, Us_pre, Rs_pre);

            total_dv = h/2*(dv_left_left + dv_right_left);
            gradR(network_left, j) = gradR(network_left, j) + total_dv;

            dv_left_right  = dv_W(F.A(i, :), F.K(i, :), j, 2, h, 0, phi_p, R1, R2, logR, Us_pre, Rs_pre);
            dv_right_right = dv_W(F.A(i, :), F.K(i, :), j, 2, h, 1, phi_p, R1, R2, logR, Us_pre, Rs_pre);

            total_dv = h/2*(dv_left_right + dv_right_right);
            gradR(network_right, j) = gradR(network_right, j) + total_dv;
        end
    end 
end
    
function [gradphi,gradR] = enforce_boudaries(F, gradphi, gradR)
    for i = 1:size(F.fixed_nodes, 2)
        index = F.network_index(F.fixed_nodes(i));
        gradphi(index, :) = zeros(1, 3);
        gradR(index, :) = zeros(1, 3);
    end
end

function L = disc_laplacian(F)
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

%Not working currently
function L = disc_curvature_laplacian(h, m, R)
    %Define discrete Laplacian
    %ASSUME grid is uniform
    %"Unpack" gradient into one column.
    L = zeros(3*m);
    

    %============= Compute block diagonal matrix action
    %============= corresponding with cross product of velocity of curve.
    %============= And an approximate action of the derivative.
    S = zeros(3*m);
    S_p = zeros(3*m);
    for i = 0:(m-1)
        if i == 0 || i== m-1
            angle_axis = quat2axisangle(R(i+1, :));
            angle_axis_p = angle_axis;
        else
            %angle_axis = 1/2*(quat2axisangle(R(i, :)) + quat2axisangle(R(i+1, :)));
            angle_axis = quat2axisangle(R(i+1, :));
            angle_axis_p = quat2axisangle(R(i+1, :)) - quat2axisangle(R(i, :));
        end

        %Index into 3x3 block diagonal
        S((3*i+1):(3*i+3), (3*i+1):(3*i+3)) = angle_axis(1)*axis2matrix(angle_axis(2:4));
        S_p((3*i+1):(3*i+3), (3*i+1):(3*i+3)) = angle_axis_p(1)*axis2matrix(angle_axis_p(2:4));
    end

    %===========Add cross product action
    L = L + S;

    %============Add approx 1st deriv
    off_diag = 1/2*ones(3*m-3,1);
    L = L + 2*diag(off_diag, 3);
    L = L - 2*diag(off_diag, -3);

    %Fix start and end... (Do this?)
    %L(1:3, 1:3) = L(1:3, 1:3) - 2*eye(3);
    %L(1:3, 4:6) = L(1:3, 4:6) + eye(3);

    %L((end-2):end, (end-2):end) = L((end-2):end, (end-2):end) + 2*eye(3);
    %L((end-2):end, (end-5):(end-3)) = L((end-2):end, (end-5):(end-3)) - eye(3);

    %============Act by cross-product
    L = S*L;

    %============Add action by derivative of curve
    L = L + S_p;

    %============Add approx 2nd deriv (OBS. bad at start and end...)
    main_diag = -2*ones(3*m,1);
    off_diag = 1*ones(3*m-3,1);
    L = L + diag(main_diag);
    L = L + diag(off_diag, 3);
    L = L + diag(off_diag, -3);

    %============Multiply by scale factor 1/h^2 and negate.
    L = -1*L/(h*h);
end



