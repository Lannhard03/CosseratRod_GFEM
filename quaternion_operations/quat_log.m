function axisangle = quat_log(R1, R2)
    %Compute log and convert to axis angle.
    %Note. Matlab still copies R1 since that is the default behaviour :(
    R1(2:4) = -R1(2:4);
    axisangle = quat_mul(R1, R2);
    %R1(2:4) = -R1(2:4); Would be neccesary if matlab didn't copy..
    
    axisangle = quat2axisangle(axisangle);
end