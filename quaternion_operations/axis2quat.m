function R = axis2quat(R_axis)
    R = zeros(1, 4);
    angle = norm(R_axis);
    axis = R_axis / angle;
    R(1) = cos(angle / 2);
    R(2:4) = sin(angle / 2) * axis.';
end