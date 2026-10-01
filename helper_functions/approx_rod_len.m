function len = approx_rod_len(phi)
    m = size(phi);
    m = m(1, 1);
    len = 0;
    for i = 1:(m-1)
        len = len + norm(phi(i, :) - phi(i+1, :)); 
    end
end