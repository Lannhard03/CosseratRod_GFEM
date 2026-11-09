function axisangle = quat_log(R1, R2)
    %Compute log and convert to axis angle.
    %Each row represent one quaternion.

    %Conjugate
    R1(:, 2:4) = -R1(:, 2:4);
    axisangle = quat_mul(R1, R2);

    %Swap sign to convention.
    axisangle(axisangle(:, 1) < 0) = -axisangle(axisangle(:, 1) < 0);
    
    
    axisangle = quat2axisangle(axisangle);
    %axisangle(1) = axisangle(1)/sqrt(2);
end