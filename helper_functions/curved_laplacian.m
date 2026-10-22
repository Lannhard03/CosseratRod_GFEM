function L = curved_laplacian(F)
    %Need an element for each direction.
    n_net = F.num_network_index;
    L = zeros(3 * n_net);
    
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
    
        n_left  = F.network_index(left);
        n_right = F.network_index(right);
    
        h = F.G(i);
    
        R1 = F.R(left, :);
        R2 = F.R(right, :);
    
        logR = quat_log(R1, R2);
        L_i = logR(1) * logR(2:4); %Convert to tangent vector 
    
        L_norm = norm(L_i);
        if L_norm == 0
            %Limit case of the below recovers naive laplacian.
            for j = 0:2
                if n_left == 0 && n_right == 0
                    continue;
                elseif n_left == 0 %In this case network_right != 0
                    L(n_right + n_net*j, n_right + n_net*j) = L(n_right + n_net*j, n_right + n_net*j)  - 1/h;
                elseif n_right == 0 %in this case network_left != 0
                    L(n_left  + n_net*j, n_left  + n_net*j) = L(n_left  + n_net*j, n_left  + n_net*j)  - 1/h;
                else
                    L(n_right + n_net*j, n_right + n_net*j) = L(n_right + n_net*j, n_right + n_net*j)  - 1/h;
                    L(n_left  + n_net*j, n_left  + n_net*j) = L(n_left  + n_net*j, n_left  + n_net*j)  - 1/h;
    
                    L(n_left  + n_net*j, n_right + n_net*j) = L(n_left  + n_net*j, n_right + n_net*j)  + 1/h;
                    L(n_right + n_net*j, n_left  + n_net*j) = L(n_right + n_net*j, n_left  + n_net*j)  + 1/h;
                end
            end
        else 
    
            K = eye(3) - ((L_i/h).'*(L_i/h))/(L_norm/h)^2;
            S = (L_i/h).'*(L_i/h)/(L_norm/h)^2;
    
            b = L_norm/(2*sqrt(2));
            s = sin(b);
            c = cos(b);
    
            main = -K/(h*s^2) * b/2 * (b + s*c) - S/h;
            side =  K/(h*s^2) * b/2 * (s + b*c) + S/h;
    
    
            for j = 0:2
                for l = 0:2
                    if n_left == 0 && n_right == 0
                        continue;
                    elseif n_left == 0 %In this case network_right != 0
                        L(n_right + n_net*j, n_right + n_net*j) = L(n_right + n_net*j, n_right + n_net*j)  + main(j+1, l+1);
                    elseif n_right == 0 %in this case network_left != 0
                        L(n_left  + n_net*j, n_left  + n_net*j) = L(n_left  + n_net*j, n_left  + n_net*j)  + main(j+1, l+1);
                    else
                        L(n_right + n_net*j, n_right + n_net*j) = L(n_right + n_net*j, n_right + n_net*j)  + main(j+1, l+1);
                        L(n_left  + n_net*j, n_left  + n_net*j) = L(n_left  + n_net*j, n_left  + n_net*j)  + main(j+1, l+1);
    
                        L(n_left  + n_net*j, n_right + n_net*j) = L(n_left  + n_net*j, n_right + n_net*j)  + side(j+1, l+1);
                        L(n_right + n_net*j, n_left  + n_net*j) = L(n_right + n_net*j, n_left  + n_net*j)  + side(j+1, l+1);
                    end
                end
            end
        end
    end
end