%Test that network indicies are computed correctly.

F = FiberNetwork;

F.m = 10;

F.fixed_nodes = [5, 9];

F.connected_nodes = {[1, 4, 3], [2, 7]};

F = F.apply_connections();

% Network indicies should be
correct = [1, 2, 1, 1, 0, 3, 2, 4, 0, 5]

F.network_index

if F.network_index == correct
    fprintf("Network indexing seems to be working\n")
else
    fprintf("Network indexing seems to be wrong\n")
end

