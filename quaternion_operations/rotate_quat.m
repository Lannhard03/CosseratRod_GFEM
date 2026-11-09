function phi_rot = rotate_quat(phi, R)
    %Each row of phi is assumed to be a 3d vector.
    %Each row of R is assumed to be a quaternion.

    c = R(:, 1);
    phi_rot = phi;

    %Don't rotate if c = 1
    phi_c = phi(c ~= 1, :);
    if isempty(phi_c) %Can we do this branchlessly???
        return;
    end

    s = sqrt(1-c(c ~= 1).^2);  %We don't care about the sign of sine here.
    l = R(c ~= 1, 2:4) ./ s;   %as the sign cancels here!

    %Rotation formula for quaternions
    phi_rot(c ~= 1, :) = phi_c + 2*c(c ~= 1).*s.*cross(l, phi_c, 2) + 2*s.^2.*cross(l, cross(l, phi_c, 2), 2);
end