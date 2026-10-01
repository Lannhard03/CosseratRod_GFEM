function F = network_solver(F, params)
    rods = F.rods;
    n = size(rods);
    n = n(1,2);

    m_tot = 0;
    J0 = 0; 

    for i = 1:n
        m_tot = m_tot + rods(i).m;
        J0 = J0 + J(rods(i).G, rods(i).A, rods(i).K, rods(i).phi, rods(i).R, rods(i).Us_pre, rods(i).Rs_pre);    
    end
    
    gradphi = zeros(m_tot, 3);
    gradR = zeros(m_tot, 3);

    h = rods(1).G(2) - rods(1).G(1);
    discrete_laplacian = disc_laplacian(F, n, h, m_tot);

    failed_backtracks = 0;
    %===== Outer loop =====
    for p = 1:params.iter_max
        fprintf("Iteration %d\n", p);

        %Compute total gradient
        curr_index = 1;
        for i = 1:n
            l = curr_index;
            r = curr_index + rods(i).m - 1;
            [gradphi(l:r, :), gradR(l:r, :)] = grad_J(rods(i).G, rods(i).A, rods(i).K, rods(i).phi, rods(i).R, rods(i).Us_pre, rods(i).Rs_pre);
            curr_index = r + 1;
        end

        %Pre-conditioning R with Laplacian gives worse results,
        %what is the correct way to pre-condition it?
        %discrete_quat_laplacian = disc_quat_laplacian(h, m_tot, rods(1).R);
        %gradR_col = reshape(gradR.', [], 1);
        %gradR_col = linsolve(discrete_quat_laplacian, gradR_col);
        %gradR = reshape(gradR_col, 3, []).';
       
        [gradphi,gradR] = enforce_connections(F,gradphi,gradR);

        %Pre-condition with discrete laplacian
        gradphi = linsolve(discrete_laplacian, gradphi);

        [gradphi,gradR] = enforce_boudaries(F,gradphi,gradR);
    
        %Compute length of gradient
        gradphi_len = 0;
        gradR_len = 0;
        for i = 1:m_tot
            gradphi_len = gradphi_len + norm(gradphi(i, :));
            %WHY 1/2 here??
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

            phi1 = zeros(m_tot, 3);
            R1 = zeros(m_tot, 4);

            %Take step in gradient direction and update energy.
            curr_index = 1;
            J1 = 0; 
            for i = 1:n
                a = curr_index;
                b = curr_index + rods(i).m - 1;
                [phi1(a:b, :), R1(a:b, :)] = step(rods(i).phi, rods(i).R, -gradphi(a:b, :), -gradR(a:b, :), steplength);

                J1 = J1 + J(rods(i).G, rods(i).A, rods(i).K, phi1(a:b, :), R1(a:b, :), rods(i).Us_pre, rods(i).Rs_pre);    
                curr_index = b + 1;
            end

            %Update if new energy is lower
            if (J0 - J1) >= params.sigma*steplength*(gradphi_len + gradR_len)

                curr_index = 1;
                for i = 1:n
                    a = curr_index;
                    b = curr_index + rods(i).m - 1;

                    rods(i).phi = phi1(a:b, :);
                    rods(i).R = R1(a:b, :);

                    curr_index = b + 1;
                end

                break;
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

        F.rods = rods;
    end
end

function [phi1, R1] = step(phi, R, gradphi, gradR, steplength)
    m = size(phi);
    m = m(1, 1);
    phi1 = phi + steplength * gradphi;
    R1 = zeros(m, 4); %allocates new mem, bad...
    gradR = steplength*gradR;
    for l = 1:m
        angle = norm(gradR(l, :));
        if angle ~= 0
            axis = gradR(l, :) / angle;
            R1(l, :) = quat_axisangle_mul(R(l, :), axis, angle);
        else
            R1(l, :) = R(l, :);
        end
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

function J = J(G, A, K, phi, R, Us_pre, Rs_pre)
    J = 0;
    m = size(G);
    m = m(1, 2);
    
    %Trapezoid quadrature.
    for i = 1:(m-1)
        phi1 = phi(i, :);
        phi2 = phi(i+1, :);
    
        R1 = R(i, :);
        R2 = R(i+1, :);
    
        h = G(i+1) - G(i);
    
        Rs = Rs_approx(h, R1, R2);
    
        %Left endpoint
        Us_pre_i = Us_pre(i, :);
        Rs_pre_i = Rs_pre(i, :);
    
        Us = Us_approx(h, phi1, phi2, R1);
    
        extension_term = dot(A, (Us - Us_pre_i).^2);
        curvature_term = 4*dot(K, (Rs - Rs_pre_i).^2);
    
        J = J + h/2*(extension_term + curvature_term);
    
        %Right endpoint
        Us_pre_i = Us_pre(i+1, :);
        Rs_pre_i = Rs_pre(i+1, :);
    
        Us = Us_approx(h, phi1, phi2, R2);
    
        extension_term = dot(A, (Us - Us_pre_i).^2);
        curvature_term = 4*dot(K, (Rs - Rs_pre_i).^2);

        J = J + h/2*(extension_term + curvature_term);
    end
end

function [gradphi, gradR] = grad_J(G, A, K, phi, R, Us_pre, Rs_pre)
    %Assemble each grid interval,
    %Do the simple thing and integrate with a trapezoid
    %quadrature rule, so that we only have to evaluate 
    %integral at nodal points.
    m = size(G);
    m = m(1, 2);
    
    gradphi = zeros(m, 3);
    gradR = zeros(m, 3);
    
    %Compute for interior nodes
    for i = 2:(m-1)
        phi1 = phi(i-1, :);
        phi2 = phi(i, :);
        phi3 = phi(i+1, :);
    
        R1 = R(i-1, :);
        R2 = R(i, :);
        R3 = R(i+1, :);
    
        h1 = G(i) - G(i-1);
        h2 = G(i+1) - G(i);

        phi1_p = (phi2 - phi1) / h1;
        phi2_p = (phi3 - phi2) / h2;
    
        logR1 = quat_log(R1, R2);
        logR2 = quat_log(R2, R3);

        Us_pre_1 = Us_pre(i-1, :);
        Rs_pre_1 = Rs_pre(i-1, :);
        Us_pre_2 = Us_pre(i, :);
        Rs_pre_2 = Rs_pre(i, :);
        Us_pre_3 = Us_pre(i+1, :);
        Rs_pre_3 = Rs_pre(i+1, :);
    
        for j = 1:3
            dw_left1 = dw_W(A, j, 2, h1, 0, phi1_p, R1, logR1, Us_pre_1);
            dw_right1 = dw_W(A, j, 2, h1, 1, phi1_p, R1, logR1, Us_pre_2);
    
            dw_left2 = dw_W(A, j, 1, h2, 0, phi2_p, R2, logR2, Us_pre_2);
            dw_right2 = dw_W(A, j, 1, h2, 1, phi2_p, R2, logR2, Us_pre_3);
    
            dw_total = h1/2*(dw_left1 + dw_right1) + h2/2*(dw_left2 + dw_right2);
            gradphi(i, j) = gradphi(i, j) + dw_total; 
    

            dv_left1 = dv_W(A, K, j, 2, h1, 0, phi1_p, R1, R2, logR1, Us_pre_1, Rs_pre_1);
            dv_right1 = dv_W(A, K, j, 2, h1, 1, phi1_p, R1, R2, logR1, Us_pre_2, Rs_pre_2);
    
            dv_left2 = dv_W(A, K, j, 1, h2, 0, phi2_p, R2, R3, logR2, Us_pre_2, Rs_pre_2);
            dv_right2 = dv_W(A, K, j, 1, h2, 1, phi2_p, R2, R3, logR2, Us_pre_3, Rs_pre_3);
    
            total_dv = h1/2*(dv_left1 + dv_right1) + h2/2*(dv_left2 + dv_right2);
            gradR(i, j) = gradR(i, j) + total_dv;
        end
    end

    %Left endpoint
    phi2 = phi(1, :);
    phi3 = phi(2, :);

    R2 = R(1, :);
    R3 = R(2, :);

    h2 = G(2) - G(1);

    phi2_p = (phi3 - phi2) / h2;

    logR2 = quat_log(R2, R3);

    Us_pre_2 = Us_pre(1, :);
    Rs_pre_2 = Rs_pre(1, :);
    Us_pre_3 = Us_pre(2, :);
    Rs_pre_3 = Rs_pre(2, :);

    for j = 1:3
        dw_left2 = dw_W(A, j, 1, h2, 0, phi2_p, R2, logR2, Us_pre_2);
        dw_right2 = dw_W(A, j, 1, h2, 1, phi2_p, R2, logR2, Us_pre_3);

        dw_total = h2/2*(dw_left2 + dw_right2);
        gradphi(1, j) = gradphi(1, j) + dw_total; 

        dv_left2 = dv_W(A, K, j, 1, h2, 0, phi2_p, R2, R3, logR2, Us_pre_2, Rs_pre_2);
        dv_right2 = dv_W(A, K, j, 1, h2, 1, phi2_p, R2, R3, logR2, Us_pre_3, Rs_pre_3);

        total_dv = h2/2*(dv_left2 + dv_right2);
        gradR(1, j) = gradR(1, j) + total_dv;
    end

    %Right endpoint
    phi1 = phi(m-1, :);
    phi2 = phi(m, :);

    R1 = R(m-1, :);
    R2 = R(m, :);

    h1 = G(m) - G(m-1);

    phi1_p = (phi2 - phi1) / h1;

    logR1 = quat_log(R1, R2);

    Us_pre_1 = Us_pre(m-1, :);
    Rs_pre_1 = Rs_pre(m-1, :);
    Us_pre_2 = Us_pre(m, :);
    Rs_pre_2 = Rs_pre(m, :);

    for j = 1:3
        dw_left1 = dw_W(A, j, 2, h1, 0, phi1_p, R1, logR1, Us_pre_1);
        dw_right1 = dw_W(A, j, 2, h1, 1, phi1_p, R1, logR1, Us_pre_2);

        dw_total = h1/2*(dw_left1 + dw_right1);
        gradphi(m, j) = gradphi(m, j) + dw_total; 

        dv_left1 = dv_W(A, K, j, 2, h1, 0, phi1_p, R1, R2, logR1, Us_pre_1, Rs_pre_1);
        dv_right1 = dv_W(A, K, j, 2, h1, 1, phi1_p, R1, R2, logR1, Us_pre_2, Rs_pre_2);

        total_dv = h1/2*(dv_left1 + dv_right1);
        gradR(m, j) = gradR(m, j) + total_dv;
    end
    

end

function [gradphi,gradR] = enforce_connections(F,gradphi,gradR)
    num_connections = size(F.connected_nodes);
    num_connections = num_connections(1);
    
    comb_gradphi = zeros(1, 3);
    comb_gradR = zeros(1, 3);
    for i = 1:num_connections
        connected = F.connected_nodes(i, :);
        num_connected = size(connected);
        num_connected = num_connected(2);
        %Sum up gradients
        for j = 1:num_connected
            comb_gradphi = comb_gradphi + gradphi(connected(j), :);
            comb_gradR = comb_gradR + gradR(connected(j), :);
        end
        %Assign combined gradient
        for j = 1:num_connected
            gradphi(connected(j), :) = comb_gradphi;
            gradR(connected(j), :) = comb_gradR;
        end
    end
end
    
function [gradphi,gradR] = enforce_boudaries(F,gradphi,gradR)
    num_fixed = size(F.fixed_nodes);
    num_fixed = num_fixed(2);
    for i = 1:num_fixed
        gradphi(F.fixed_nodes(i), :) = zeros(1, 3);
        gradR(F.fixed_nodes(i), :) = zeros(1, 3);
    end
end

function L = disc_laplacian(F, n, h, m_tot)
    L = zeros(m_tot);

    num_connections = size(F.connected_nodes);
    num_connections = num_connections(1);

    %Connections within each rod
    curr_index = 1;
    for i = 1:n
        m = F.rods(i).m;
        a = curr_index;
        b = curr_index + F.rods(i).m - 1;
        
        main_diag = -2*ones(m,1);
        off_diag = 1*ones(m-1,1);
        %Add to block corresponding with rod i
        L(a:b, a:b) = L(a:b, a:b) + diag(main_diag);
        L(a:b, a:b) = L(a:b, a:b) + diag(off_diag, 1);
        L(a:b, a:b) = L(a:b, a:b) + diag(off_diag,-1);
        
       
        curr_index = b + 1;
    end

    %Connections between rods
    for i = 1:num_connections
        connected = F.connected_nodes(i, :);
        num_connected = size(connected);
        num_connected = num_connected(2);

        for k = 1:num_connected
            for l = 1:num_connected
                if k == l
                    continue;
                end
                %L(connected(k), connected(l)) = L(connected(k), connected(l)) + 1;
            end
        end
    end

    L = -1*L/(h*h);
end

function L = disc_quat_laplacian(h, m, R)
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



