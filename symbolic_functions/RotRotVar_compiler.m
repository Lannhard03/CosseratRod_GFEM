function RotRotVar_compiler()
    %w_i = [1, 0, 1]; %(v_i - v_pre_i)^T M_2
    %E_j = [1, 0, 0]; %Test direction
    %v_i = [0, 1, 0]; %velocity

    %E_orth = E_j - dot(v_i, E_j)*v_i;
    %b_i = pi/2; %\norm(velocity)/2sqrt(2).

    syms t;
    syms b;
    syms h;
    syms E [1,3];
    syms w [1,3];
    syms v [1,3];

    E_orth = E - dot(v, E)*v;

    s_i  = sin(b * t) / sin(b) * E_orth;

    s_ip = diff(s_i, t);

    s_i_cross = cross(v, s_i);

    integrand_sign = dot(w, s_ip + E); %Derivative terms, swaps sign
    integrand = dot(w, s_i_cross); %Non-derivative term, same sign

    integral_sign = int(integrand_sign, t, 0, 1);
    integral = int(integrand, t, 0, 1);

    matlabFunction(integral_sign, 'Vars', {h, b, E, w, v}, ...
        "File", "symbolic_functions/RotRotVar_sign", "Optimize",true);
    matlabFunction(integral, 'Vars', {h, b, E, w, v}, ...
        "File", "symbolic_functions/RotRotVar", "Optimize",true);
end