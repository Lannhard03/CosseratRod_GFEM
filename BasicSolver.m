%0. Do assembly, i.e compute for each grid interval at a time...

%1. Need to: Compute gradient of energy functional (same as lifted 
%            energy functional at 0 in tangent space).

%2. Need to: Compute gradient of lifted internal energy density
%            and use quadrature to compute 1.

%3. Need to: Compute derivatives in w_ij and v_ij directions 
%            (i = 1, ..., m, j = 1, 2, 3), of the algebraic strains Us, Rs.

%4. Need to: Do usual FEM hat stuff for w_ij deriv of Us
% (what are the hat functions for values in R^3?)

%5. Need to: Compute v_ij deriv geodesic gamma[Q_i exp(v_i), Q_i+1
%exp(v_i+1)](..scaling..), to solve for v_ij deriv of Us and Rs.

%6. Derivative from 5. reduces to an integral over matrix stuff...
%(quadrature again?).

%Very inefficient as matrices are used throughout?
%Okay thats fine...

%Obs the initial values are hard coded into phi, and R
%will not be changed by alg.
function [phi, R] = solver(G, A, K, phi, R, Us_pre, Rs_pre)
    %Use backtracking method...
    m = size(G);
    m = m(1, 2);
    J0 = J(G, A, K, phi, R, Us_pre, Rs_pre);
    %Max 10 gradient descent steps,
    %exit if gradient is zero.

    alpha = 0.0001;
    beta = 0.9;
    sigma = 0.1;

    for p = 1:1000
        fprintf("Iteration %d\n", p);
        [gradphi, gradR] = grad_J(G, A, K, phi, R, Us_pre, Rs_pre);
        %Force start and end point to be fixed
        gradphi(:, 1) = zeros(3, 1);
        gradphi(:, m) = zeros(3, 1);
        gradR(:, :, 1) = zeros(3, 3);
        gradR(:, :, m) = zeros(3, 3); 
        gradphi_len = 0;
        for l = 1:m
            gradphi_len = gradphi_len + norm(gradphi(:, l))^2;
        end

        gradR_len = 0;
        for l = 1:m
            gradR_len= gradR_len + trace(gradR(:, :, l).'*gradR(:, :, l));
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
            
            [phi1, R1] = step(phi, R, -gradphi, -gradR, steplength);

            J1 = J(G, A, K, phi1, R1, Us_pre, Rs_pre);
    
            if (J0 - J1) >= sigma*steplength*(gradphi_len + gradR_len)
                phi = phi1;
                R = R1;
                if l < 500
                    sigma = min(1.1*sigma, 1);
                end
                break;
            end
        end

        J0 = J1;
        fprintf("Backtracking took %d iterations \n", l);

        if l == 1000
            fprintf("Gradient descent found no better solution iteration %d :(\n", p);
            sigma = 0.7*sigma;
            alpha = beta*alpha;
            %break;
        end

        if l == 1
            alpha = alpha/beta;
        end



        % J0 = J(G, A, K, phi, R, Us_pre, Rs_pre);
        % for l = 0:1000
        %     steplength = beta^l * alpha;
        %     [phi1, R1] = step(phi, R, 0*gradphi, -gradR, steplength);
        % 
        %     J1 = J(G, A, K, phi, R1, Us_pre, Rs_pre);
        % 
        %     if (J0 - J1) >= sigma*steplength*gradR_len
        %         R = R1;
        %         break;
        %     end
        % end
        % 
        % if l == 1000
        %     fprintf("Gradient descent found no better solution for R iteration %d :(\n", p);
        % end
    end
end

function [phi, R] = quat_solver(G, A, K, phi, R, Us_pre, Rs_pre)
    %Use backtracking method...
    m = size(G);
    m = m(1, 2);
    J0 = quat_J(G, A, K, phi, R, Us_pre, Rs_pre);
    %Max 10 gradient descent steps,
    %exit if gradient is zero.
    
    alpha = 1;
    beta = 0.9;
    sigma = 0.1;
    
    for p = 1:3000
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
                if l < 500
                    sigma = min(1.1*sigma, 1);
                end
                break;
            end
        end
    
        fprintf("Backtracking took %d iterations \n", l);
    
        if l == 1000
            fprintf("Gradient descent found no better solution iteration %d :(\n", p);
            sigma = 0.7*sigma;
            alpha = beta*alpha;
            if alpha < 0.01
                break;
            end
        else
            J0 = J1;
        end
    
        if l == 1
            alpha = alpha/beta;
        end
    
    
    
        % J0 = J(G, A, K, phi, R, Us_pre, Rs_pre);
        % for l = 0:1000
        %     steplength = beta^l * alpha;
        %     [phi1, R1] = step(phi, R, 0*gradphi, -gradR, steplength);
        % 
        %     J1 = J(G, A, K, phi, R1, Us_pre, Rs_pre);
        % 
        %     if (J0 - J1) >= sigma*steplength*gradR_len
        %         R = R1;
        %         break;
        %     end
        % end
        % 
        % if l == 1000
        %     fprintf("Gradient descent found no better solution for R iteration %d :(\n", p);
        % end
    end
end


function [phi1, R1] = step(phi, R, gradphi, gradR, steplength)
    m = size(phi);
    m = m(1, 2);
    phi1 = phi + steplength * gradphi;
    R1 = zeros(3,3,m); %allocates new mem, bad...
    for l = 1:m
        R1(:, :, l) = R(:, :, l)*expm(steplength * gradR(:, :, l));
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

% Keep everything as quaternions, no conversion to angle...
% Renormalize after multiplication?
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
        dir(j) = -1;
    else
        dir(j) = 1;
    end
    R1_cong = R1;
    %Change sign of imaginary part
    R1_cong(2:4) = -1*R1_cong(2:4);
    aux = quat2axisangle(quat_mul(R2, R1_cong));
    ei = rotate_quat(dir, R1);
    l = aux(2:4);
    a = aux(1);
    %if the angle is zero the derivative is zero
    if a == 0
        deriv = 0*l;
        return;
    end
    c = cos(a);
    s = sin(a);

    deriv = 1/h*(a*c/s * ei + dot(l, ei)*(1 - a*c/s)*l - abs(a)*cross(l, ei));
    %When going to the "tangent space" we must scale by 2/pi !
    %deriv = 2*aux(1)/(pi*h) * aux(2:4);
end

%Compute dv_Us for one element at
%the left endpoint if k = 1, and the right if k = 2.
function deriv = quat_dv_Us(j, k, h, phi_p, R1, R2)
    dir = zeros(1,4);
    dir(j+1) = -1;
    %Maybe uneccesary to create this aux variable?

    if k == 1
        R = R1;
    elseif k == 2
        R = R2;
    end
    
     R(2:4) = -R(2:4); %conjugate

    aux1 = rotate_quat(phi_p, R);
    quat = zeros(1,4);
    quat(2:4) = aux1;


    quat = 1/(2) * (quat_mul(dir,quat) - quat_mul(quat, dir));

    deriv = quat(2:4);
end


function deriv = quat_dw_Us(j, k, h, s, R1, R2)
    dir = zeros(1,3);
    if k == 1
        dir(j) = -1/h;
    elseif k == 2
        dir(j) = 1/h;
    end
    
    %Compute log and convert to axis angle.
    R1_cong = R1;
    R1_cong(2:4) = -R1_cong(2:4);
    aux = quat_mul(R1_cong, R2);
    aux = quat2axisangle(aux);
    
    %scale rotation and reverse direction
    aux(1) = -2*s*aux(1);
    dir = rotate_axisangle(dir, aux(1), aux(2:4));
    
    deriv = rotate_quat(dir, R1_cong);
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

function approx = quat_dv_W(A, K, j, k, h, s, phi1, phi2, R1, R2, Us_pre, Rs_pre)
    %Contributions from extension strain
    R1_cong = R1;
    R1_cong(2:4) = -R1_cong(2:4);
    quat = quat_mul(R2, R1_cong);
    angle_axis = quat2axisangle(quat);
    phi_p = 1/h*(phi2 - phi1);

    %Rotate first by exp(s*log(R1'*R2)), then by R1.
    aux = rotate_axisangle(phi_p, -2*s*angle_axis(1), angle_axis(2:4));
    Us = rotate_quat(aux, R1_cong);
    %Are we compting Us wrong as R*phi_p instead of as R'*phi_p??

    factor1 = A.*(Us - Us_pre);
    

    %Not compute for general s in [0,1]...
    factor2 = quat_dv_Us(j, k, h, phi_p, R1, R2); 
    approx = dot(factor1, factor2);
    
    %Contributions from curvature strain
    Rs_coords = 2*quat_Rs_approx(h, R1, R2);
    
    scaled_diff= K.*(Rs_coords - Rs_pre);
    
    dv = quat_dv_Rs(j, k, h, R1, R2);
    
    approx = approx + dot(scaled_diff, dv);
end

function approx = quat_dw_W(A, j, k, h, s, phi1, phi2, R1, R2, Us_pre)
    %Contributions from extension strain
    R1_cong = R1;
    R1_cong(2:4) = -R1_cong(2:4);
    quat = quat_mul(R2, R1_cong);
    angle_axis = quat2axisangle(quat);
    phi_p = 1/h*(phi2 - phi1);
    
    %Rotate first by exp(s*log(R1'*R2)), then by R1.
    aux = rotate_axisangle(phi_p, -2*s*angle_axis(1), angle_axis(2:4));
    Us = rotate_quat(aux, R1_cong);
    
    factor1 = A.*(Us - Us_pre);
    
    factor2 = quat_dw_Us(j, k, h, s, R1, R2);
    
    approx = dot(factor1, factor2);
end

function strain = Rs_approx(h, R1, R2)
    Rs_mat = 1/h*logm(R1.' * R2);
    strain = [Rs_mat(3,2), Rs_mat(1,3), Rs_mat(2,1)];
end

function strain = Us_approx(h, phi1, phi2, R)
    phi_p = 1/h*(phi2 - phi1);
    strain = R.' * phi_p;
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
        curvature_term = dot(K, (Rs - Rs_pre_i).^2);
    
        J = J + h/2*(extension_term + curvature_term);
    
        %Right endpoint
        Us_pre_i = Us_pre(i+1, :);
        Rs_pre_i = Rs_pre(i+1, :);
    
        Us = quat_Us_approx(h, phi1, phi2, R2);
    
        extension_term = dot(A, (Us - Us_pre_i).^2);
        curvature_term = dot(K, (Rs - Rs_pre_i).^2);
    
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
    
    %We zero the gradient for the exterior nodes so ignore them
    %The gradient is computed node wise.
    for i = 2:(m-1)
        phi1 = phi(i-1, :);
        phi2 = phi(i, :);
        phi3 = phi(i+1, :);
    
        R1 = R(i-1, :);
        R2 = R(i, :);
        R3 = R(i+1, :);
    
        h1 = G(i) - G(i-1);
        h2 = G(i+1) - G(i);
    
        Us_pre_1 = Us_pre(i-1, :);
        Rs_pre_1 = Rs_pre(i-1, :);
        Us_pre_2 = Us_pre(i, :);
        Rs_pre_2 = Rs_pre(i, :);
        Us_pre_3 = Us_pre(i+1, :);
        Rs_pre_3 = Rs_pre(i+1, :);
    
        for j = 1:3
            dw_left1 = quat_dw_W(A, j, 2, h1, 0, phi1, phi2, R1, R2, Us_pre_1);
            dw_right1 = quat_dw_W(A, j, 2, h1, 1, phi1, phi2, R1, R2, Us_pre_2);
    
            dw_left2 = quat_dw_W(A, j, 1, h2, 0, phi2, phi3, R2, R3, Us_pre_2);
            dw_right2 = quat_dw_W(A, j, 1, h2, 1, phi2, phi3, R2, R3, Us_pre_3);
    
            dw_total = h1/2*(dw_left1 + dw_right1) + h2/2*(dw_left2 + dw_right2);
            gradphi(i, j) = gradphi(i, j) + dw_total; 
    
    
            dv_left1 = quat_dv_W(A, K, j, 2, h1, 0, phi1, phi2, R1, R2, Us_pre_1, Rs_pre_1);
            dv_right1 = quat_dv_W(A, K, j, 2, h1, 1, phi1, phi2, R1, R2, Us_pre_2, Rs_pre_2);
    
            dv_left2 = quat_dv_W(A, K, j, 1, h2, 0, phi2, phi3, R2, R3, Us_pre_2, Rs_pre_2);
            dv_right2 = quat_dv_W(A, K, j, 1, h2, 1, phi2, phi3, R2, R3, Us_pre_3, Rs_pre_3);
    
            total_dv = h1/2*(dv_left1 + dv_right1) + h2/2*(dv_left2 + dv_right2);
            gradR(i, j) = gradR(i, j) + total_dv;
        end
    end

end

function J = J(G, A, K, phi, R, Us_pre, Rs_pre)
    J = 0;
    m = size(G);
    m = m(1, 2);

    %Trapezoid quadrature.
    for i = 1:(m-1)
        phi1 = phi(:, i);
        phi2 = phi(:, i+1);

        R1 = R(:, :, i);
        R2 = R(:, :, i+1);

        h = G(i+1) - G(i);

        Rs_coords = Rs_approx(h, R1, R2);

        %Left endpoint
        Us_pre_i = Us_pre(:, i);
        Rs_pre_i = Rs_pre(:, :, i);

        Us = Us_approx(h, phi1, phi2, R1);
        Rs_pre_i_coords = [Rs_pre_i(3,2), Rs_pre_i(1,3), Rs_pre_i(2,1)];

        extension_term = dot(A, (Us - Us_pre_i).^2);
        curvature_term = dot(K, (Rs_coords - Rs_pre_i_coords).^2);

        J = J + h/2*(extension_term + curvature_term);

        %Right endpoint
        Us_pre_i = Us_pre(:, i+1);
        Rs_pre_i = Rs_pre(:, :, i+1);

        Us = Us_approx(h, phi1, phi2, R2);
        Rs_pre_i_coords = [Rs_pre_i(3,2), Rs_pre_i(1,3), Rs_pre_i(2,1)];

        extension_term = dot(A, (Us - Us_pre_i).^2);
        curvature_term = dot(K, (Rs_coords - Rs_pre_i_coords).^2);

        J = J + h/2*(extension_term + curvature_term);
    end
end


%G is vector of start and end points of the grid intervals
%i.e 0 = l_1, l_2, ..., l_m = 1.
%A, K material parameters,
%phi, R are algebraic nodal values,
%i.e a list of m vectors in R^3, and a list of rotations resp.
%Us_pre, Rs_pre, are algebraic pre-strains for the extension
%and curvature strains
function [gradphi, gradR] = grad_J(G, A, K, phi, R, Us_pre, Rs_pre)
    %Assemble each grid interval,
    %Do the simple thing and integrate with a trapezoid
    %quadrature rule, so that we only have to evaluate 
    %integral at nodal points.
    m = size(G);
    m = m(1, 2);
       
    gradphi = zeros(3, m);
    gradR = zeros(3, 3, m);

    %We zero the gradient for the exterior nodes so ignore them
    %The gradient is computed node wise.
    for i = 2:(m-1)
        phi1 = phi(:, i-1);
        phi2 = phi(:, i);
        phi3 = phi(:, i+1);

        R1 = R(:, :, i-1);
        R2 = R(:, :, i);
        R3 = R(:, :, i+1);

        h1 = G(i) - G(i-1);
        h2 = G(i+1) - G(i);

        Us_pre_1 = Us_pre(:, i-1);
        Rs_pre_1 = Rs_pre(:, :, i-1);
        Us_pre_2 = Us_pre(:, i);
        Rs_pre_2 = Rs_pre(:, :, i);
        Us_pre_3 = Us_pre(:, i+1);
        Rs_pre_3 = Rs_pre(:, :, i+1);

        for j = 1:3
            dw_left1 = dw_W(A, j, 2, h1, 0, phi1, phi2, R1, R2, Us_pre_1);
            dw_right1 = dw_W(A, j, 2, h1, 1, phi1, phi2, R1, R2, Us_pre_2);

            dw_left2 = dw_W(A, j, 1, h2, 0, phi2, phi3, R2, R3, Us_pre_2);
            dw_right2 = dw_W(A, j, 1, h2, 1, phi2, phi3, R2, R3, Us_pre_3);

            %Terrible can cancellations happen here...
            dw_total = h1/2*(dw_left1 + dw_right1) + h2/2*(dw_left2 + dw_right2);
            gradphi(j, i) = gradphi(j, i) + dw_total; 


            dv_left1 = dv_W(A, K, j, 2, h1, 0, phi1, phi2, R1, R2, Us_pre_1, Rs_pre_1);
            dv_right1 = dv_W(A, K, j, 2, h1, 1, phi1, phi2, R1, R2, Us_pre_2, Rs_pre_2);


            dv_left2 = dv_W(A, K, j, 1, h2, 0, phi2, phi3, R2, R3, Us_pre_2, Rs_pre_2);
            dv_right2 = dv_W(A, K, j, 1, h2, 1, phi2, phi3, R2, R3, Us_pre_3, Rs_pre_3);

             %Terrible can cancellations happen here...
            total_dv = h1/2*(dv_left1 + dv_right1) + h2/2*(dv_left2 + dv_right2);
            gradR(:, :, i) = gradR(:, :, i) + total_dv*skew_basis(j);
        end
    end


end

function approx = dw_W(A, j, k, h, s, phi1, phi2, R1, R2, Us_pre)
    %Contributions from extension strain
    R = R1*expm(s * logm(R1.' * R2));
    Us = Us_approx(h, phi1, phi2, R);
    
    factor1 = A.'.*(Us - Us_pre);

    factor2 = dw_Us(j, k, h, s, R1, R2);
    
    approx = dot(factor1, factor2);
end

function approx = dv_W(A, K, j, k, h, s, phi1, phi2, R1, R2, Us_pre, Rs_pre)
    %Contributions from extension strain
    R = R1*expm(s * logm(R1.' * R2));
    phi_p = 1/h*(phi2 - phi1);
    Us = R.' * phi_p;
    factor1 = A.'.*(Us - Us_pre);
    
    %This ought to be wrong, derivative on the i'th node
    %should depend on the positions of the (i-1)'th and (i+1)'th node.
    factor2 = dv_Us(j, k, h, s, phi_p, R1, R2); 
    approx = dot(factor1, factor2);
    
    %Contributions from curvature strain
    Rs_coords = Rs_approx(h, R1, R2);

    Rs_pre_coords = [Rs_pre(3,2), Rs_pre(1,3), Rs_pre(2,1)];
    scaled_diff_coords = K.*(Rs_coords - Rs_pre_coords);
    
    %This ought to be wrong, derivative on the i'th node
    %should depend on the positions of the (i-1)'th and (i+1)'th node.
    dv = dv_Rs(j, k, h, R1, R2);
    dv_coords = [dv(3,2), dv(1,3), dv(2,1)];
    
    approx = approx + dot(scaled_diff_coords, dv_coords);
end

%Approximate derivative of Us in the w_jk direction,
%at a value s.
%k = 1, 2 (left or right grid point of grid interval)
%j = 1, 2, 3 (direction in R^3)
%h grid size
%R1, R2 in SO(3),
%s point along rod affinely mapped to the standard element [0,1].
function deriv = dw_Us(j, k, h, s, R1, R2)
    if k == 1
        hat_deriv = -1/h;
    elseif k == 2
        hat_deriv = 1/h;
    end

    R = R1*expm(s * logm(R1.' * R2));
    
    %The j'th row from each column.
    %Test with transponant here..
    deriv = hat_deriv * R(j, :); 
end

%Approximate derivative of Us in the w_jk direction,
%at a value s.
%k = 1, 2 (left or right grid point of grid interval)
%j = 1, 2, 3 (direction in R^3)
%h grid size
%R1, R2 in SO(3),
%s point along rod affinely mapped to the standard element [0,1].
function approx = dv_Us(j, k, h, s, phi_p, R1, R2)
    delta = h/10000; %Step for finite difference... good, bad??
    E_j = skew_basis(j);
    step = expm(delta*E_j);
    if k == 1
        step_inside = step.';
        step_outside = step;
    elseif k == 2
        step_inside = step.';
        step_outside = step;
        s = 1-s;

        R = R1;
        R1 = R2;
        R2 = R; 
    end

    M1 = R1*expm(s * logm(R1.' * R2));
    M2 = R1*step_outside*expm(s * logm(R1.' * step_inside * R2));
    deriv = 1/(delta) * (M2 - M1);

    approx = deriv.' * phi_p; %d/dvjk(R^T) * phi' ??
end


%Approximate derivative of Rs in the v_jk direction.
%use a stupid right rectangle rule...
%k = 1, 2 (left or right grid point of grid interval)
%j = 1, 2, 3 (direction in so(3))
%h grid size
%R1, R2 in SO(3),
function approx = dv_Rs(j, k, h, R1, R2)
    delta = h/1000; %Step for finite difference... good, bad??
    E_j = skew_basis(j);
    if k == 1
        step = expm(delta*E_j).';
    elseif k == 2
        step = expm(delta*E_j);
    end

    M1 = logm(R1.' * R2);
    M2 = logm(R1.' * step * R2);
    approx = 1/(h*delta) * (M2 - M1);
end

%return the k'th basis matrix for the skew-symmetric matrices
function E_k = skew_basis(k)
    E_k = zeros(3,3);
    if k == 1
        E_k(2, 3) = -1;
        E_k(3, 2) = 1;
    end
    if k == 2
        E_k(1, 3) = 1;
        E_k(3, 1) = -1;
    end
    if k == 3
        E_k(1, 2) = -1;
        E_k(2, 1) = 1;
    end
end


function draw_rod(phi, R)
    figure;
    xlim([-0.3 1.5])
    ylim([-0.9 0.9])
    zlim([-0.9 0.9])
    view(3);
    hold on
    
    m = size(phi);
    m = m(1, 2);
    
    for i = 1:(m)
        base = phi(:, i);
        R1 = R(:, :, i) / 10;
        quiver3(base(1), base(2), base(3), R1(1, 1), R1(2, 1), R1(3, 1), 'red', 'LineWidth', 2)
        quiver3(base(1), base(2), base(3), R1(1, 2), R1(2, 2), R1(3, 2), 'blue', 'LineWidth', 2)
        quiver3(base(1), base(2), base(3), R1(1, 3), R1(2, 3), R1(3, 3), 'green', 'LineWidth', 2)
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


%Initialize solution as the direct
%geodesic interpolation between a fixed start and end point.
function [phi, R] = interpolation_start(m, phi1, phi2, R1, R2)
    G = (0:(m-1))/(m-1);
    %Set up rod:
    I = eye(3, 3);
    %Weird web code to create 3d matrix...
    R = repmat(I,[1,1,m]);
    
    phi = phi1 .* (1-G) + phi2 .* G;

    Rs = logm(R1.' * R2);
    for l = 1:m
        R(:,:,l) = R1 * expm(G(l)*Rs);
    end
end

function len = approx_rod_len(phi)
    m = size(phi);
    m = m(1, 2);
    len = 0;
    for i = 1:(m-1)
       len = len + norm(phi(:, i) - phi(:, i+1)); 
    end
end


m = 32;
G = (0:(m-1))/(m-1);
%Set up rod:
phi1 = [0, 0, 0].';
phi2 = [0.5, 0, 0].';





a = pi/2; 
R1_q = [cos(a/2), 0, 0, sin(a/2)];
R1 = [  cos(a), -sin(a), 0;
        sin(a),  cos(a), 0;
             0,       0, 1;
     ];

b = pi/2; 
R2 = [ cos(b), 0, sin(b);
            0, 1,      0;
      -sin(b), 0, cos(b);
     ];
R2_q = [cos(b/2), 0, sin(b/2), 0];

%phi2 = R1*[1, 0, 0].';

[phi, R] = interpolation_start(m, phi1, phi2, R1, R2);
init_len = approx_rod_len(phi);
phi = zeros(m, 3);
phi(:, 1) = 0.5*(0:(m-1))/(m-1); %Constant velocity in x-direction

%draw_rod(phi, R);

A = [1963, 755, 755];
K = [0.94, 1.23, 1.23];

%Pre strains at grid points
%representing a straight rod of length 1.
Us_pre = zeros(3, m);
Us_pre(1, :) = ones([1, m]); %Constant velocity in x-direction
Rs_pre = zeros(3, 3, m); %No rotation to begin with.
%Rs_pre(2,3,:) = ones([1, m]);
%Rs_pre(3,2,:) = -ones([1, m]); %Pre-strain rotation in the direction of the rod...

Rs_pre_quat = zeros(m, 3);
Us_pre_quat = Us_pre.';

h = 1;
k = 2;
j = 1;


phi1 = [0, 3, 0];
phi2 = [2,1,1];


interp = linspace(0, 1, m);
R_q = zeros(m, 4);
for l = 1:m
    R_q(l, :) = interp(l)*R2_q + (1-interp(l))*R1_q;
    R_q(l, :) = R_q(l, :) / norm(R_q(l, :));
end


%angles = linspace(a/2, b/2, m);
%R_q = zeros(m, 4);
%R_q(:, 1) = cos(angles);
%R_q(:, 4) = sin(angles);

%Compare rotation matrix and quaternion implementations.

%====== TEST dv_Rs!!! ======  Seems to work!
%deriv = quat_dv_Rs(j, k, h, R1_q, R2_q)
%approx = dv_Rs(j, k, h, R1, R2)


%======  TEST dv_Us!!! ======  Seems to work!
% phi_p = [2,0,0];
% s = k-1;
% approx = dv_Us(j, k, h, s, phi_p.', R1, R2)
% deriv = quat_dv_Us(j, k, h, phi_p, R1_q, R2_q)


%======  TEST dw_Us!!! ======  Seems to work!
% s = k-1;
% approx = dw_Us(j, k, h, s, R1, R2)
% deriv = quat_dw_Us(j, k, h, s, R1_q, R2_q)


%======  TEST dv_W!!! ====== Seems to work!
% s = k-1;
% deriv = quat_dv_W(A, K, j, k, h, s, phi1, phi2, R1_q, R2_q, Us_pre_quat(1, :), Rs_pre_quat(1, :))
% approx = dv_W(A, K, j, k, h, s, phi1.', phi2.', R1, R2, Us_pre(:, 1), Rs_pre(:, :,1))


%======  TEST dw_W!!! ====== Seems to work!
% s = k-1;
% approx = dw_W(A, j, k, h, s, phi1.', phi2.', R1, R2, Us_pre(:, 1))
% deriv = quat_dw_W(A, j, k, h, s, phi1, phi2, R1_q, R2_q, Us_pre_quat(1, :))


%======  TEST energy!!! ====== Seems to work!
% energy = quat_J(G, A, K, phi, R_q, Us_pre_quat, Rs_pre_quat)
% energy_approx = J(G, A, K, phi.', R, Us_pre, Rs_pre)

%======  TEST gradient of energy!!! ====== Seems to work!
% [gradphi_q, gradR_q] = quat_grad_J(G, A, K, phi, R_q, Us_pre_quat, Rs_pre_quat);
% [gradphi, gradR] = grad_J(G, A, K, phi.', R, Us_pre, Rs_pre);
% gradphi_q
% gradphi


%======  TEST step!!! ====== Seems to work!
% E_j = 5*pi*skew_basis(1);
% [phi3, R3] = step(phi1.', R1, [1,0,0].', E_j, 0.1)
% 
% gradR_q = 5*pi*[1,0,0];
% [phi3_quat, R3_quat] = quat_step(phi1, R1_q, [1,0,0], gradR_q, 0.1);
%phi3_quat
%quat2axisangle(R3_quat)


%======  TEST step with gradient!! ======
% [gradphi_q, gradR_q] = quat_grad_J(G, A, K, phi, R_q, Us_pre_quat, Rs_pre_quat);
% [phi3_quat, R3_quat] = quat_step(phi, R_q, gradphi_q, gradR_q, 0.0001);
% 
% [gradphi, gradR] = grad_J(G, A, K, phi.', R, Us_pre, Rs_pre);
% [phi3, R3] = step(phi.', R, gradphi, gradR, 0.0001);
% phi3_quat
% phi3

%[phi, R] = solver(G, A, K, phi, R, Us_pre, Rs_pre);
%end_len = approx_rod_len(phi);
%fprintf("Rod started at %d lenght, and ended at %d, should have length 1\n", init_len, end_len);


quat_draw_rod(phi, R_q);

[phi, R] = quat_solver(G, A, K, phi, R_q, Us_pre_quat, Rs_pre_quat);

quat_draw_rod(phi, R);
% 
% draw_rod(phi.', R);
% 
% [phi, R] = solver(G, A, K, phi.', R, Us_pre, Rs_pre);
% 
% draw_rod(phi, R);