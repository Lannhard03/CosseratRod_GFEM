function mat = axis2matrix(l)
    %convert axis to skew-symmetric matrix
    %cross product by l, becomes matrix multiplication.
    mat = zeros(3);
    
    mat(1, 2) = -l(3);
    mat(2, 1) =  l(3);
    
    mat(1, 3) =  l(2);
    mat(3, 1) = -l(2);
    
    mat(2, 3) = -l(1);
    mat(3, 2) =  l(1);
end