function phi_rot = rotate_axisangle(phi, angle, axis)
    c = cos(angle);
    s = sin(angle);
    phi_rot = c*phi +  (1-c)*dot(axis,phi)*axis + s*cross(axis, phi);
end