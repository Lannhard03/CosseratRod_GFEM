classdef FiberNetwork
    properties
        rods
        %OBS!!! The connected and fixed nodes need to be disjoint!!!
        fixed_nodes {mustBeNumeric}
        connected_nodes {mustBeNumeric}
    end
end
