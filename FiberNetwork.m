classdef FiberNetwork
    properties
        %rods
        %OBS!!! The connected and fixed nodes need to be disjoint!!!
        fixed_nodes {mustBeNumeric}

        %Each row must be sorted in increasing order
        connected_nodes

        m {mustBeNumeric} %Total number of nodes
        phi {mustBeNumeric}
        R {mustBeNumeric}
        edges {mustBeNumeric}

        %Map indicies to the correct ones in the network (where nodes are
        %conected).
        num_network_index {mustBeNumeric}
        network_index {mustBeNumeric}

        G {mustBeNumeric} % h/reference length for each edge.
        %Material param for each edge
        A {mustBeNumeric}
        K {mustBeNumeric} 

        %Pre strain for each edge
        U_pre {mustBeNumeric}
        R_pre {mustBeNumeric}

        %Dependent quantities
        theta {mustBeNumeric}
        n     {mustBeNumeric}
        phi_p {mustBeNumeric}
        Us   {mustBeNumeric}
        Rs   {mustBeNumeric}
    end

    methods
        %Compute network indexing
        function F = apply_connections(F)
            num_groups = size(F.connected_nodes, 2);
            total_connected = 0;
            for i = 1:num_groups
                total_connected = total_connected + size(F.connected_nodes{i}, 2);
            end

            num_fixed = size(F.fixed_nodes, 2);
            
            %Remove the number of indicies not needed.
            F.num_network_index = F.m - (total_connected - num_groups) - num_fixed;

            F.network_index = 1:F.m;
            %Compensate indicies after a connected node index
            for i = 1:num_groups
                group = F.connected_nodes{i};
                for j = 2:size(group, 2)
                    F.network_index(group(j):end) = F.network_index(group(j):end) - 1;
                end
            end

            %Compensate indicies after a fixed node index
            for i = 1:num_fixed
                F.network_index(F.fixed_nodes(i):end) = F.network_index(F.fixed_nodes(i):end) - 1;
                F.network_index(F.fixed_nodes(i)) = 0;
            end

            %Set correct index for each connected group
            for i = 1:num_groups
                group = F.connected_nodes{i};
                for j = 2:size(group, 2)
                    F.network_index(group(j)) =  F.network_index(group(1));
                end
            end
        end
       
        %Update dependent quantities of the rod
        function F = update_network(F)
            left = F.edges(:, 1);
            right = F.edges(:, 2);


            phi_left  = F.phi(left, :);
            phi_right = F.phi(right, :);
            F.phi_p = (phi_right - phi_left) ./ F.G;


            R_left  = F.R(left, :);
            R_right = F.R(right, :);

            logR = quat_log(R_left, R_right);
            F.theta = logR(:, 1);
            F.n = logR(:, 2:4);
            F.Rs = logR(:, 1) .* logR(:, 2:4) ./ F.G; %Convert to tangent vector

            R_left(:, 2:4) = -R_left(:, 2:4);
            F.Us = rotate_quat(F.phi_p, R_left);      
        end
    end
end
