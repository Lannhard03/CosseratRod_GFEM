function rot = quat_mul(R1, R2)
    %Each row is treated as one quaternion
    s1 = R1(:, 1);
    s2 = R2(:, 1);
    n1 = R1(:, 2:4);
    n2 = R2(:, 2:4);
    
    rot = R1;
    rot(:, 1) = s1.*s2 - dot(n1, n2, 2);
    rot(:, 2:4) = s2.*n1 + s1.*n2 + cross(n1, n2, 2);
end