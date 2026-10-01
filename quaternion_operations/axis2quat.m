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