function angle_axis = quat2axisangle(R)
    %Each row is assumed to be one quaternion.

    %The first entry in the output is the angle of the rotation
    %and the three last are the rotation axis.
    %OBS. quaternions form a double cover,
    %     inverting rotation angle and axis gives the same rotation.

    angle_axis = zeros(size(R));
    angle_axis(:, 1) = acos(R(:, 1));
    s = sin(angle_axis(:, 1));
    if isempty(s(s ~= 0)) % Can we do this branchlessly??
        return; 
    end

    angle_axis(s ~= 0, 2:4) = R(s ~= 0, 2:4)./s(s ~= 0);
    %The remaining ones are only zero.
end