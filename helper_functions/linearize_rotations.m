function r = linearize_rotations(R)
    %R is a mx4 matrix where each row is a quaternion.
    %r is a mx3 matrix where each row is the linearization of R

    temp = quat2axisangle(R);
    r = 2*temp(:, 1) .* temp(:, 2:4);
end