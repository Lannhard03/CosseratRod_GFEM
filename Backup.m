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

    %Max 10 gradient descent steps,
    %exit if gradient is zero.
    for p = 1:10
        [gradphi, gradR] = grad_J(G, A, K, phi, R, Us_pre, Rs_pre);
        %Force start and end point to be fixed
        gradphi(:, 1) = zeros(3, 1);
        gradphi(:, m) = zeros(3, 1);
        gradR(:, :, 1) = zeros(3, 3);
        gradR(:, :, m) = zeros(3, 3);
        grad_len = 0;
        for l = 1:m
            grad_len = grad_len + norm(gradphi(:, l))^2 + norm(gradR(:, :, l))^2;
        end
    
        alpha = 1;
        beta = 0.5;
        sigma = 0.5;
        
        %should be a while loop here...
        J0 = J(G, A, K, phi, R, Us_pre, Rs_pre);
        for l = 0:100
            steplength = beta^l * alpha;
            [phi1, R1] = step(phi, R, -gradphi, -gradR, steplength);

            J1 = J(G, A, K, phi1, R1, Us_pre, Rs_pre);
    
            if (J0 - J1) >= sigma*steplength*grad_len
                phi = phi1;
                R = R1;
                break;
            end
        end

        if l == 100
            fprintf("Gradient descent found no better solution :(\n");
        end

        if grad_len == 0
           break;
        end
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


function J = J(G, A, K, phi, R, Us_pre, Rs_pre)
    J = 0;
    m = size(G);
    m = m(1, 2);
    
    for i = 1:(m-1)
        phi1 = phi(:, i);
        phi2 = phi(:, i+1);

        R1 = R(:, :, i);
        R2 = R(:, :, i+1);

        h = G(i+1) - G(i);

        Us_pre_i = Us_pre(:, i);
        Rs_pre_i = Rs_pre(:, :, i);
        
        phi_p = -1/h*phi1 + 1/h*phi2;
        Us = R1.' * phi_p;
        Rs = 1/h*logm(R1.' * R2);

        Rs_coords = [Rs(3,2), Rs(1,3), Rs(2,1)];
        Rs_pre_i_coords = [Rs_pre_i(3,2), Rs_pre_i(1,3), Rs_pre_i(2,1)];

        extension_term = dot(A, (Us - Us_pre_i).^2);
        curvature_term = dot(K, (Rs_coords - Rs_pre_i_coords).^2);
        
        %We use trapezoid quadrature so endpoints
        %have half the weight.
        if i == 1 || i == m
            factor = 1/4;
        else
            factor = 1/2;
        end

        J = J + factor*(extension_term + curvature_term);
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

    for i = 1:(m-1)
        phi1 = phi(:, i);
        phi2 = phi(:, i+1);

        R1 = R(:, :, i);
        R2 = R(:, :, i+1);

        h = G(i+1) - G(i);

        Us_pre_i = Us_pre(:, i);
        Rs_pre_i = Rs_pre(:, :, i);

        for k = 1:2
            for j = 1:3
                dw = dw_W(A, j, k, h, 0, phi1, phi2, R1, R2, Us_pre_i);
                dv = dv_W(A, K, j, k, h, 0, phi1, phi2, R1, R2, Us_pre_i, Rs_pre_i);

                gradphi(j, i-1+k) = gradphi(j, i-1+k) + dw;

                gradR_skew = zeros(3,3);
                for q = 1:3
                    E_q = skew_basis(q);
                    gradR_skew = gradR_skew + dv(q)*E_q;
                end
                

                gradR(:, :, i-1+k) = gradR(:, :, i-1+k) + gradR_skew;
            end
        end
    end


end


function approx = dw_W(A, j, k, h, s, phi1, phi2, R1, R2, Us_pre)
    %Contributions from extension strain
    phi_p = -1/h*phi1 + 1/h*phi2;
    R = R1*expm(s * logm(R1.' * R2));
    Us = R.' * phi_p;
    factor1 = A.'.*(Us - Us_pre);

    %This ought to be wrong, derivative on the i'th node
    %should depend on the positions of the (i-1)'th and (i+1)'th node.
    factor2 = dw_Us(j, k, h, s, R1, R2);
    
    approx = dot(factor1, factor2);
end

function approx = dv_W(A, K, j, k, h, s, phi1, phi2, R1, R2, Us_pre, Rs_pre)
    %Contributions from extension strain
    phi_p = -1/h*phi1 + 1/h*phi2;
    R = R1*expm(s * logm(R1.' * R2));
    Us = R.' * phi_p;
    factor1 = A.'.*(Us - Us_pre);
    
    %This ought to be wrong, derivative on the i'th node
    %should depend on the positions of the (i-1)'th and (i+1)'th node.
    factor2 = dv_Us(j, k, h, s, phi_p, R1, R2); 
    approx = dot(factor1, factor2);
    
    %Contributions from curvature strain
    Rs = 1/h * logm(R1.' * R2);
    factor1 = K.*(Rs - Rs_pre);
    
    %This ought to be wrong, derivative on the i'th node
    %should depend on the positions of the (i-1)'th and (i+1)'th node.
    factor2 = dv_Rs(j, k, h, R1, R2);
    
    approx = approx + dot(factor1, factor2);
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
    delta = h/1000; %Step for finite difference... good, bad??
    E_j = skew_basis(j);
    if k == 1
        step = expm(delta*E_j).';
    elseif k == 2
        step = expm(delta*E_j);
    end

    M1 = R1*expm(s * logm(R1.' * R2));
    M2 = R1*expm(s * logm(R1.' * step * R2));
    deriv = 1/delta * (M2 - M1);

    approx = deriv.' * phi_p; %d/dvjk(R^T) * phi' ??
end


%Approximate derivative of Rs in the v_jk direction.
%use a stupid right rectangle rule...
%k = 1, 2 (left or right grid point of grid interval)
%j = 1, 2, 3 (direction in so(3))
%h grid size
%R1, R2 in SO(3),
%There is some symmetry for s and 1-s... unused
%Better to just use a simple f.d rule?
function approx = dv_Rs(j, k, h, R1, R2)
    approx = zeros(3,3);
    if k == 1
        sign = -1;
    else
        sign = 1;
    end

    E_j = skew_basis(j);

    n = 10; %10 quadrature points, too few? enough?
    for i = 1:n
        s = i/n;
        factor_1 = (s*R2 + (1-s)*R1)^-1;
        factor_2 = (s*R1 + (1-s)*R2)^-1;
        term = 1/n * factor_1 * E_j * factor_2.';
        approx = approx + term;
    end

    approx = 1/h * sign * approx;
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
    m = m(1, 2)
    
    for i = 1:(m)
        base = phi(:, i);
        R1 = R(:, :, i) / 10;
        for j = 1:3
            quiver3(base(1), base(2), base(3), R1(j, 1), R1(j, 2), R1(j, 3), 'blue', 'LineWidth', 2)
        end
    end
end


m = 10;
G = (0:(m-1))/(m-1);
%Set up rod:
I = eye(3, 3);
%Weird web code to create 3d matrix...
R = repmat(I,[1,1,m]);

phi = zeros(3, m);
x_coord = (0:(m-1))/(m-1);
phi(1, :) = x_coord;

phi(:, m) = [0.5, 0.3, 0];

a = pi/2;
R(:, :, m) = [ cos(a),  sin(a), 0;
              -sin(a),  cos(a), 0;
                    0,       0, 1;
             ];    

draw_rod(phi, R);

A = 0.1*[1963, 755, 755];
K = [1.23, 1.23, 0.94];

Us_pre = zeros(3, m);
Us_pre(1, :) = ones([1, m]); %Constant velocity in x-direction
Rs_pre = zeros(3, 3, m); %No rotation to begin with.


[phi, R] = solver(G, A, K, phi, R, Us_pre, Rs_pre);

draw_rod(phi, R);
J(G, A, K, phi, R, Us_pre, Rs_pre)