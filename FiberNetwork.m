classdef FiberNetwork
    properties
        %rods
        %OBS!!! The connected and fixed nodes need to be disjoint!!!
        fixed_nodes {mustBeNumeric}

        %Each row must be sorted in increasing order
        connected_nodes {mustBeNumeric}

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
        Us_pre {mustBeNumeric}
        Rs_pre {mustBeNumeric}
    end

    methods
        function F = apply_connections(F)
            %Wrong for more connected groups!!!
            total_connections = size(F.connected_nodes, 2);
            connected_groups = size(F.connected_nodes, 1);
            F.num_network_index = F.m - total_connections + connected_groups;

            F.network_index = 1:F.m;
            for i = 1:size(F.connected_nodes, 1)
                connected = F.connected_nodes(i, :);
                for j = 2:size(connected, 2)
                    F.network_index(connected(j):end) = F.network_index(connected(j):end) - 1;
                    F.network_index(connected(j)) =  F.network_index(connected(1));
                end
            end
        end
       

    end
end
