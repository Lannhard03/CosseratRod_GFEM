function F_l = normal2linearNetwork(F)
F_l = LinearFiberNetwork;

F_l.m = F.m;
F_l.G = F.G;

F_l.edges = F.edges;

F_l.A = F.A;
F_l.K = F.K;

F_l.fixed_nodes = F.fixed_nodes;
F_l.connected_nodes = F.connected_nodes;

F_l = F_l.apply_connections();

%Initial displacement is zero.
F_l.u = zeros(F.m, 3);

%Initial positions taken from phi.
F_l.start = F.phi;

F_l.r = linearize_rotations(F.R);

%Edge direction n is the pre-strain.
F_l.n = F.U_pre;
end