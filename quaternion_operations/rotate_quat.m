function phi_rot = rotate_quat(phi, R)
    c = R(1);
    if c == 1
        phi_rot = phi;
        return;
    end
    s = sqrt(1-c^2); %We dont care about the sign of sine,
    l = R(2:4) / s;  %as the sign cancels here!
    
    %Rotation formula for quaternions
    phi_rot = phi + 2*c*s*cross(l, phi) + 2*s^2*cross(l, cross(l, phi));
end