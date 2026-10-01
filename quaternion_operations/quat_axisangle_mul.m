function rot = quat_axisangle_mul(R, axis, angle)
    s1 = R(1);
    s2 = cos(angle/2);
    n1 = R(2:4);
    n2 = sin(angle/2)*axis;
    
    rot = zeros(1, 4);
    rot(1) = s1*s2 - dot(n1, n2);
    rot(2:4) = s2*n1 + s1*n2 + cross(n1, n2);
end