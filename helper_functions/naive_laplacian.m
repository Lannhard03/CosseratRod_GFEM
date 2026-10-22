function L = naive_laplacian(F)
    L = zeros(F.num_network_index);
    
    for i = 1:size(F.edges, 1)
        left = F.edges(i, 1);
        right = F.edges(i, 2);
    
        network_left  = F.network_index(left);
        network_right = F.network_index(right);
    
        h = F.G(i);
    
        %Cases: 
        % Both left and right fixed
        % Left fixed
        % Right fixed
        % No fixed
        if network_left == 0 && network_right == 0
            continue;
        elseif network_left == 0 %In this case network_right != 0
            L(network_right, network_right)= L(network_right, network_right) - 1/h;
        elseif network_right == 0 %in this case network_left != 0
            L(network_left, network_left)  = L(network_left, network_left)   - 1/h;
        else
            L(network_left, network_left)  = L(network_left, network_left)   - 1/h;
            L(network_right, network_right)= L(network_right, network_right) - 1/h;
    
            L(network_left, network_right) = L(network_left, network_right)  + 1/h;
            L(network_right, network_left) = L(network_right, network_left)  + 1/h;
        end
    end
end