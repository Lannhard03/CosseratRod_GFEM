function R_mat = quat2mat(R)    
    m = size(R);
    m = m(1, 1);
    R_mat = zeros(3, 3, m); 
    for i = 1:(m)
        R_mat(:, 1, i) = rotate_quat([1,0,0], R(i, :));
        R_mat(:, 2, i) = rotate_quat([0,1,0], R(i, :));
        R_mat(:, 3, i) = rotate_quat([0,0,1], R(i, :));
    end
end