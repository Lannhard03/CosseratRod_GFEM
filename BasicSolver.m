%TODO: Maybe combine all rods into combined memory phi, R and only keep
%track of start and endpoints for each rod (needed to compute correct
%gradients and energy.

function F = network_solver(F)
    beta = 0.9;
    alpha = 0.01;
    sigma = 1e-6;

    rods = F.rods;
    n = size(rods);
    n = n(1,2);

    m_tot = 0;
    J0 = 0; 

    for i = 1:n
        m_tot = m_tot + rods(i).m;
        J0 = J0 + quat_J(rods(i).G, rods(i).A, rods(i).K, rods(i).phi, rods(i).R, rods(i).Us_pre, rods(i).Rs_pre);    
    end
    
    gradphi = zeros(m_tot, 3);
    gradR = zeros(m_tot, 3);
    
    for p = 1:1000
        fprintf("Iteration %d\n", p);

        m_curr = 1;
        for i = 1:n
            l = m_curr;
            r = m_curr + rods(i).m - 1;
            [gradphi(l:r, :), gradR(l:r, :)] = quat_grad_J(rods(i).G, rods(i).A, rods(i).K, rods(i).phi, rods(i).R, rods(i).Us_pre, rods(i).Rs_pre);
            m_curr = r + 1;
        end
    
        gradR = gradR/2; %Weird scale factor!?!
        
        %Enforce fixed points
        num_fixed = size(F.fixed_nodes);
        num_fixed = num_fixed(2);
        for i = 1:num_fixed
            gradphi(F.fixed_nodes(i), :) = zeros(1, 3);
            gradR(F.fixed_nodes(i), :) = zeros(1, 3);
        end

        %Enforce connected points
        %I.e make sure the gradients for the points are the same
        %Are they just addable or do we need to transform them to the same
        %frame??
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

        gradphi_len = 0;
        gradR_len = 0;
        for i = 1:m_tot
            gradphi_len = gradphi_len + norm(gradphi(l, :))^2;
            gradR_len= gradR_len + 1/2 * norm(gradR(l, :));
        end
    
        fprintf("gradient lengths are %d and %d\n", gradphi_len, gradR_len);
    
        fprintf("Energy is %d\n", J0);
        if gradR_len + gradphi_len < 10^(-15)
            break;
        end
    
        for l = 0:1000
            steplength = beta^l * alpha;

            phi1 = zeros(m_tot, 3);
            R1 = zeros(m_tot, 4);

            m_curr = 1;
            J1 = 0; 
            for i = 1:n
                a = m_curr;
                b = m_curr + rods(i).m - 1;
                [phi1(a:b, :), R1(a:b, :)] = quat_step(rods(i).phi, rods(i).R, -gradphi(a:b, :), -gradR(a:b, :), steplength);

                J1 = J1 + quat_J(rods(i).G, rods(i).A, rods(i).K, phi1(a:b, :), R1(a:b, :), rods(i).Us_pre, rods(i).Rs_pre);    
                m_curr = b + 1;
            end

            if (J0 - J1) >= sigma*steplength*(gradphi_len + gradR_len)

                m_curr = 1;
                for i = 1:n
                    a = m_curr;
                    b = m_curr + rods(i).m - 1;

                    rods(i).phi = phi1(a:b, :);
                    rods(i).R = R1(a:b, :);

                    m_curr = b + 1;
                end

                break;
            end
        end
    
        fprintf("Backtracking took %d iterations \n", l);
    
        if l == 1000
            fprintf("Gradient descent found no better solution iteration %d :(\n", p);
    
        else
            J0 = J1;
        end
    
        if l <= 300
            %sigma = 0.7*sigma;
            %alpha = alpha/beta;
        end
    
        F.rods = rods;
    end
end


function [phi, R] = quat_solver(G, A, K, phi, R, Us_pre, Rs_pre)
    %Use backtracking method...
    m = size(G);
    m = m(1, 2);
    J0 = quat_J(G, A, K, phi, R, Us_pre, Rs_pre);    
    
    beta = 0.5;
    alpha = 0.1;
    sigma = 1e-8;
    
    for p = 1:1000
        fprintf("Iteration %d\n", p);
        [gradphi, gradR] = quat_grad_J(G, A, K, phi, R, Us_pre, Rs_pre);
        gradR = gradR/2; %Weird scale factor!?!
        %Force start and end point to be fixed
        gradphi(1, :) = zeros(1, 3);
        gradphi(m, :) = zeros(1, 3);
        gradR(1, :) = zeros(1, 3);
        gradR(m, :) = zeros(1, 3); 
        gradphi_len = 0;
        for l = 1:m
            gradphi_len = gradphi_len + norm(gradphi(l, :))^2;
        end
    
        gradR_len = 0;
        for l = 1:m
            gradR_len= gradR_len + 1/2 * norm(gradR(l, :));
        end
    
    
        fprintf("gradient lengths are %d and %d\n", gradphi_len, gradR_len);
    
        fprintf("Energy is %d\n", J0);
        if gradR_len + gradphi_len < 10^(-15)
            break;
        end
    
    
    
        %should be a while loop here...
        %separate steps in phi and R since
        %stuff doesn't work...
    
        for l = 0:1000
            steplength = beta^l * alpha;

            [phi1, R1] = quat_step(phi, R, -gradphi, -gradR, steplength);

            J1 = quat_J(G, A, K, phi1, R1, Us_pre, Rs_pre);

            if (J0 - J1) >= sigma*steplength*(gradphi_len + gradR_len)
                phi = phi1;
                R = R1;
                if l > 300
                    %sigma = min(1.1*sigma, 1);
                    %alpha = beta*alpha;
                end
                break;
            end
        end
            
    
        fprintf("Backtracking took %d iterations \n", l);
    
        if l == 1000
            fprintf("Gradient descent found no better solution iteration %d :(\n", p);

        else
            J0 = J1;
        end
    
        if l <= 300
            %sigma = 0.7*sigma;
            %alpha = alpha/beta;
        end
    end
end

function [phi1, R1] = quat_step(phi, R, gradphi, gradR, steplength)
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

function phi_rot = rotate_quat(phi, R)
    c = R(1);
    if c == 1
        phi_rot = phi;
        return;
    end
    s = sqrt(1-c^2); %We dont care about the sign of sine,
    l = R(2:4) / s;  %as the sign cancels here!

    phi_rot = phi + 2*c*s*cross(l, phi) + 2*s^2*cross(l, cross(l, phi));
end

function phi_rot = rotate_axisangle(phi, angle, axis)
    c = cos(angle);
    s = sin(angle);
    phi_rot = c*phi +  (1-c)*dot(axis,phi)*axis + s*cross(axis, phi);
end

function rot = quat_mul(R1, R2)
    s1 = R1(1);
    s2 = R2(1);
    n1 = R1(2:4);
    n2 = R2(2:4);

    rot = zeros(1, 4);
    rot(1) = s1*s2 - dot(n1, n2);
    rot(2:4) = s2*n1 + s1*n2 + cross(n1, n2);
end

function rot = quat_axisangle_mul(R, axis, angle)
    s1 = R(1);
    s2 = cos(angle/2);
    n1 = R(2:4);
    n2 = sin(angle/2)*axis;

    rot = zeros(1, 4);
    rot(1) = s1*s2 - dot(n1, n2);
    rot(2:4) = s2*n1 + s1*n2 + cross(n1, n2);
end

%Convert quaternion to a 4d vector,
%the first entry is the angle of the rotation
%and the three last are the rotation axis.
%OBS. quaternions form a double cover,
%     inverting rotation angle and axis gives the same rotation.
function angle_axis = quat2axisangle(R)
    angle_axis = zeros(1,4);
    angle_axis(1) = acos(R(1));
    s = sin(angle_axis(1));
    if s ~= 0
        angle_axis(2:4) = R(2:4)/s;
    else
        angle_axis(2:4) = R(2:4); %This will be zero.
    end
end


%k = 1,2 (left or right end point)
function deriv = quat_dv_Rs(j, k, h, R1, R2)
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
    %When going to the "tangent space" we must scale by 2/pi !
    %deriv = 2*aux(1)/(pi*h) * aux(2:4);
end

%Compute dv_Us for one element at
%the left endpoint if k = 1, and the right if k = 2.
function deriv = quat_dv_Us(j, k, phi_p, R1, R2)
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

function axisangle = quat_log(R1, R2)
    %Compute log and convert to axis angle.
    %Note. Matlab still copies R1 since that is the default behaviour :(
    R1(2:4) = -R1(2:4);
    axisangle = quat_mul(R1, R2);
    %R1(2:4) = -R1(2:4); Would be neccesary if matlab didn't copy..

    axisangle = quat2axisangle(axisangle);
end

function deriv = quat_dw_Us(j, k, h, s, R1, logR)
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

function strain = quat_Rs_approx(h, R1, R2)
    R1_cong = R1;
    R1_cong(2:4) = -R1_cong(2:4);
    strain_quat = quat_mul(R1_cong, R2);
    aux = quat2axisangle(strain_quat);
    strain = 1/h * aux(1) * aux(2:4);
end

function strain = quat_Us_approx(h, phi1, phi2, R)
    phi_p = 1/h*(phi2 - phi1);
    R_cong = R;
    R_cong(2:4) = -R_cong(2:4);
    
    strain = rotate_quat(phi_p, R_cong);
end

function approx = quat_dv_W(A, K, j, k, h, s, phi_p, R1, R2, logR, Us_pre, Rs_pre)
    %Take conjugate, then swap back
    R1(2:4) = -R1(2:4);
    dir = rotate_quat(phi_p, R1);
    R1(2:4) = -R1(2:4);
    
    %scale rotation by 2s and reverse direction
    Us = rotate_axisangle(dir, -2*s*logR(1), logR(2:4));

    factor1 = A.*(Us - Us_pre);

    %Not compute for general s in [0,1]...
    factor2 = quat_dv_Us(j, k, phi_p, R1, R2); 
    approx = dot(factor1, factor2);
    
    %Contributions from curvature strain
    Rs_coords = quat_Rs_approx(h, R1, R2);
    
    scaled_diff= K.*(Rs_coords - Rs_pre);
    
    dv = quat_dv_Rs(j, k, h, R1, R2);
    
    approx = approx + dot(scaled_diff, dv);
end

function approx = quat_dw_W(A, j, k, h, s, phi_p, R1, logR, Us_pre)
    %Take conjugate, then swap back
    R1(2:4) = -R1(2:4);
    aux = rotate_quat(phi_p, R1);
    R1(2:4) = -R1(2:4);

    Us = rotate_axisangle(aux, -2*s*logR(1), logR(2:4));

    factor1 = A.*(Us - Us_pre);
    
    factor2 = quat_dw_Us(j, k, h, s, R1, logR);
    
    approx = dot(factor1, factor2);
end

function J = quat_J(G, A, K, phi, R, Us_pre, Rs_pre)
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
    
        Rs = quat_Rs_approx(h, R1, R2);
    
        %Left endpoint
        Us_pre_i = Us_pre(i, :);
        Rs_pre_i = Rs_pre(i, :);
    
        Us = quat_Us_approx(h, phi1, phi2, R1);
    
        extension_term = dot(A, (Us - Us_pre_i).^2);
        curvature_term = 4*dot(K, (Rs - Rs_pre_i).^2);
    
        J = J + h/2*(extension_term + curvature_term);
    
        %Right endpoint
        Us_pre_i = Us_pre(i+1, :);
        Rs_pre_i = Rs_pre(i+1, :);
    
        Us = quat_Us_approx(h, phi1, phi2, R2);
    
        extension_term = dot(A, (Us - Us_pre_i).^2);
        curvature_term = 4*dot(K, (Rs - Rs_pre_i).^2);

        J = J + h/2*(extension_term + curvature_term);
    end
end

function [gradphi, gradR] = quat_grad_J(G, A, K, phi, R, Us_pre, Rs_pre)
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
            dw_left1 = quat_dw_W(A, j, 2, h1, 0, phi1_p, R1, logR1, Us_pre_1);
            dw_right1 = quat_dw_W(A, j, 2, h1, 1, phi1_p, R1, logR1, Us_pre_2);
    
            dw_left2 = quat_dw_W(A, j, 1, h2, 0, phi2_p, R2, logR2, Us_pre_2);
            dw_right2 = quat_dw_W(A, j, 1, h2, 1, phi2_p, R2, logR2, Us_pre_3);
    
            dw_total = h1/2*(dw_left1 + dw_right1) + h2/2*(dw_left2 + dw_right2);
            gradphi(i, j) = gradphi(i, j) + dw_total; 
    

            dv_left1 = quat_dv_W(A, K, j, 2, h1, 0, phi1_p, R1, R2, logR1, Us_pre_1, Rs_pre_1);
            dv_right1 = quat_dv_W(A, K, j, 2, h1, 1, phi1_p, R1, R2, logR1, Us_pre_2, Rs_pre_2);
    
            dv_left2 = quat_dv_W(A, K, j, 1, h2, 0, phi2_p, R2, R3, logR2, Us_pre_2, Rs_pre_2);
            dv_right2 = quat_dv_W(A, K, j, 1, h2, 1, phi2_p, R2, R3, logR2, Us_pre_3, Rs_pre_3);
    
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
        dw_left2 = quat_dw_W(A, j, 1, h2, 0, phi2_p, R2, logR2, Us_pre_2);
        dw_right2 = quat_dw_W(A, j, 1, h2, 1, phi2_p, R2, logR2, Us_pre_3);

        dw_total = h2/2*(dw_left2 + dw_right2);
        gradphi(1, j) = gradphi(1, j) + dw_total; 

        dv_left2 = quat_dv_W(A, K, j, 1, h2, 0, phi2_p, R2, R3, logR2, Us_pre_2, Rs_pre_2);
        dv_right2 = quat_dv_W(A, K, j, 1, h2, 1, phi2_p, R2, R3, logR2, Us_pre_3, Rs_pre_3);

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
        dw_left1 = quat_dw_W(A, j, 2, h1, 0, phi1_p, R1, logR1, Us_pre_1);
        dw_right1 = quat_dw_W(A, j, 2, h1, 1, phi1_p, R1, logR1, Us_pre_2);

        dw_total = h1/2*(dw_left1 + dw_right1);
        gradphi(m, j) = gradphi(m, j) + dw_total; 

        dv_left1 = quat_dv_W(A, K, j, 2, h1, 0, phi1_p, R1, R2, logR1, Us_pre_1, Rs_pre_1);
        dv_right1 = quat_dv_W(A, K, j, 2, h1, 1, phi1_p, R1, R2, logR1, Us_pre_2, Rs_pre_2);

        total_dv = h1/2*(dv_left1 + dv_right1);
        gradR(m, j) = gradR(m, j) + total_dv;
    end
    

end

function draw_network(F)
    figure;
    xlim([-0.3 1.5])
    ylim([-0.9 0.9])
    zlim([-0.9 0.9])
    view(3);
    hold on
    
    rods = F.rods;
    num_rods = size(rods);
    num_rods = num_rods(1, 2);

    for j = 1:num_rods
        rod = rods(j);
        phi = rod.phi;
        R = rod.R;

        m = size(phi);
        m = m(1, 1);

        for i = 1:m
            base = phi(i, :);
            R1 = rotate_quat([1,0,0], R(i, :)) / 10;
            R2 = rotate_quat([0,1,0], R(i, :)) / 10;
            R3 = rotate_quat([0,0,1], R(i, :)) / 10;

            quiver3(base(1), base(2), base(3), R1(1), R1(2), R1(3), 'red', 'LineWidth', 2)
            quiver3(base(1), base(2), base(3), R2(1), R2(2), R2(3), 'blue', 'LineWidth', 2)
            quiver3(base(1), base(2), base(3), R3(1), R3(2), R3(3), 'green', 'LineWidth', 2)
        end
    end
end

function quat_draw_rod(phi, R)
    figure;
    xlim([-0.3 1.5])
    ylim([-0.9 0.9])
    zlim([-0.9 0.9])
    view(3);
    hold on
    
    m = size(phi);
    m = m(1, 1);
    
    for i = 1:(m)
        base = phi(i, :);
        R1 = rotate_quat([1,0,0], R(i, :)) / 10;
        R2 = rotate_quat([0,1,0], R(i, :)) / 10;
        R3 = rotate_quat([0,0,1], R(i, :)) / 10;
    
        quiver3(base(1), base(2), base(3), R1(1), R1(2), R1(3), 'red', 'LineWidth', 2)
        quiver3(base(1), base(2), base(3), R2(1), R2(2), R2(3), 'blue', 'LineWidth', 2)
        quiver3(base(1), base(2), base(3), R3(1), R3(2), R3(3), 'green', 'LineWidth', 2)
    end
end

function R_mat = quat2mat(R)    
    m = size(R);
    m = m(1, 1);
    R_mat = zeros(3, 3, m); 
    for i = 1:(m)
        R_mat(:, 1, i) = rotate_quat([1,0,0], R(i, :));
        R_mat(:, 2, i) = rotate_quat([0,1,0], R(i, :));
        R_mat(:, 3, i) = rotate_quat([0,0,1], R(i, :));
    end
end

function len = approx_rod_len(phi)
    m = size(phi);
    m = m(1, 1);
    len = 0;
    for i = 1:(m-1)
       len = len + norm(phi(i, :) - phi(i+1, :)); 
    end
end

function R = axis2quat(R_axis)
    m = size(R_axis);
    m = m(1, 2);
    
    R = zeros(m, 4);
    for i = 1:m
        angle = norm(R_axis(:, i));
        axis = R_axis(:, i) / angle;
        R(i, 1) = cos(angle / 2);
        R(i, 2:4) = sin(angle / 2) * axis.';
    end
end


function one_rod_test_problem()
    m = 32;
    G = (0:(m-1))/(m-1);

    %Set up rod:
    phi = zeros(m, 3);
    phi(:, 1) = 0.5*(0:(m-1))/(m-1); %Constant velocity in x-direction
    
    a = pi/2; 
    b = pi/2; 

    angles = linspace(a/2, b/2, m);
    R_q = zeros(m, 4);
    R_q(:, 1) = cos(angles);
    R_q(:, 4) = sin(angles);
    
    %draw_rod(phi, R);
    
    A = [1963, 755, 755];
    K = [0.94, 1.23, 1.23];
    
    %Pre strains at grid points
    %representing a straight rod of length 1.
    Us_pre = zeros(m, 3);
    Us_pre(1, :) = ones([m, 1]); %Constant velocity in x-direction
    
    Rs_pre = zeros(m, 3);    

    %Reference solution.
    phi_ref = readmatrix("../FeniCSReferenceSolver/referenceSolPos.txt").';
    R_axis_ref = readmatrix("../FeniCSReferenceSolver/referenceSolRot.txt");
    R_ref = axis2quat(R_axis_ref);
end

function two_rod_test_problem()
    F = FiberNetwork;
    %OBS!!! The connected and fixed nodes need to be disjoint!!!

    %Start of rod 1 and end of rod 2 are fixed
    F.fixed_nodes = [1, 64]; 
    %Nodes 32 and 33 are connected: their relative positions and rotations should be preserved.
    %For more connections add more rows.
    F.connected_nodes = [32, 33]; 

    rod1 = Rod;
    rod2 = Rod;

    rod1.m = 32;
    rod1.G = linspace(0,1, rod1.m);
    rod1.A = [1963, 755, 755];
    rod1.K = [0.94, 1.23, 1.23];

    rod2.m = 32;
    rod2.G = linspace(0,1, rod2.m);
    rod2.A = [1963, 755, 755];
    rod2.K = [0.94, 1.23, 1.23];

    rod1.phi = zeros(rod1.m, 3);
    rod1.R = zeros(rod1.m, 4);

    rod1.phi(:, 1) = linspace(0, 1, 32).';
    rod1.R(:, 1) = 1;

    rod2.phi = zeros(rod2.m, 3);
    rod2.R = zeros(rod2.m, 4);

    rod2.phi(:, 1) = 1;
    %Add "bending" here!
    rod2.phi(:, 1) = linspace(1, 0.5, 32).';
    rod2.phi(:, 2) = linspace(0, 1, 32).';

    rod2.R(:, 1) = 1/sqrt(2);
    rod2.R(:, 4) = 1/sqrt(2); %Rotation 90 around z-axis

    rod1.Us_pre = zeros(rod1.m, 3);
    rod1.Us_pre(:, 1) = 1;

    rod2.Us_pre = zeros(rod2.m, 3);
    rod2.Us_pre(:, 1) = 1;

    rod1.Rs_pre = zeros(rod1.m, 3);
    rod2.Rs_pre = zeros(rod2.m, 3);

    F.rods = [rod1, rod2];

    draw_network(F);
    F = network_solver(F);
    draw_network(F)

    R1 = F.rods(1).R(32, :);
    R2 = F.rods(2).R(1, :);

    log = quat_log(R1, R2);

    phi1 = F.rods(1).phi;
    phi2 = F.rods(2).phi;
end


function F = three_rod_test_problem()
    F = FiberNetwork;
    %OBS!!! The connected and fixed nodes need to be disjoint!!!
    
    %Start of rod 1 and end of rod 2 are fixed
    F.fixed_nodes = [1, 64, 96]; 
    %Nodes 32 and 33 are connected: their relative positions and rotations should be preserved.
    %For more connections add more rows.
    F.connected_nodes = [32, 33, 65]; 
    
    rod1 = Rod;
    rod2 = Rod;
    rod3 = Rod;
    
    rod1.m = 32;
    rod1.G = linspace(0,1, rod1.m);
    rod1.A = [1963, 755, 755];
    rod1.K = [0.94, 1.23, 1.23];
    
    rod2.m = 32;
    rod2.G = linspace(0,1, rod2.m);
    rod2.A = [1963, 755, 755];
    rod2.K = [0.94, 1.23, 1.23];

    rod3.m = 32;
    rod3.G = linspace(0,1, rod3.m);
    rod3.A = [1963, 755, 755];
    rod3.K = [0.94, 1.23, 1.23];
    
    rod1.phi = zeros(rod1.m, 3);
    rod1.R = zeros(rod1.m, 4);
    
    rod1.phi(:, 1) = linspace(0, 1, 32).';
    rod1.R(:, 1) = 1;
    

    rod2.phi = zeros(rod2.m, 3);
    rod2.R = zeros(rod2.m, 4);
    
    %Add "bending" here!
    rod2.phi(:, 1) = linspace(1, 0.5, 32).';
    rod2.phi(:, 2) = linspace(0, 1, 32).';

    rod2.R(:, 1) = 1/sqrt(2);
    rod2.R(:, 4) = 1/sqrt(2); %Rotation 90 around z-axis


    rod3.phi = zeros(rod2.m, 3);
    rod3.R = zeros(rod2.m, 4);

    %Add "bending" here!
    rod3.phi(:, 1) = linspace(1, 1.5, 32).';
    rod3.phi(:, 2) = linspace(0, -1, 32).';

    rod3.R(:, 1) = 1/sqrt(2);
    rod3.R(:, 4) = -1/sqrt(2); %Rotation 90 around z-axis

    
    rod1.Us_pre = zeros(rod1.m, 3);
    rod1.Us_pre(:, 1) = 1;
    
    rod2.Us_pre = zeros(rod2.m, 3);
    rod2.Us_pre(:, 1) = 1;

    rod3.Us_pre = zeros(rod3.m, 3);
    rod3.Us_pre(:, 1) = 1;
    
    rod1.Rs_pre = zeros(rod1.m, 3);
    rod2.Rs_pre = zeros(rod2.m, 3);
    rod3.Rs_pre = zeros(rod3.m, 3);
    
    F.rods = [rod1, rod2, rod3];
end


function pres1()
    load("32gridSol_quat.mat");
    quat_draw_rod(phi, R)
    quat_draw_rod(phi_ref, R_ref);

    dist = max(vecnorm(phi - phi_ref, 2, 2))
end

function pres2()
    load("2rods_test.mat");
    draw_network(F)
end

function pres3()
    load("3rods_test.mat");
    draw_network(F)
end

%F = three_rod_test_problem();

%load("3rods_test.mat");
% draw_network(F) 
%F = network_solver(F);
%draw_network(F)
% 
% R1 = F.rods(1).R(32, :);
% R2 = F.rods(2).R(1, :);
% R3 = F.rods(3).R(1, :);
% 
% log = quat_log(R1, R3);
% 
% phi1 = F.rods(1).phi;
% phi2 = F.rods(2).phi;


pres2()
