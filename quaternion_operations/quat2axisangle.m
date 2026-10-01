function angle_axis = quat2axisangle(R)
    %Convert quaternion to a 4d vector,
    %the first entry is the angle of the rotation
    %and the three last are the rotation axis.
    %OBS. quaternions form a double cover,
    %     inverting rotation angle and axis gives the same rotation.
    angle_axis = zeros(1,4);
    angle_axis(1) = acos(R(1));
    s = sin(angle_axis(1));
    if s ~= 0
        angle_axis(2:4) = R(2:4)/s;
    else
        angle_axis(2:4) = R(2:4); %This will be zero.
    end
end