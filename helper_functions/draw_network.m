function draw_network(F)
    figure;
    xlim([-0.3 1.5])
    ylim([-0.9 0.9])
    zlim([-0.9 0.9])
    view(3);
    hold on
  
    
    m = F.m;

    for i = 1:m
        base = F.phi(i, :);
        R1 = rotate_quat([1,0,0], F.R(i, :)) / 10;
        R2 = rotate_quat([0,1,0], F.R(i, :)) / 10;
        R3 = rotate_quat([0,0,1], F.R(i, :)) / 10;

        quiver3(base(1), base(2), base(3), R1(1), R1(2), R1(3), 'red', 'LineWidth', 2)
        quiver3(base(1), base(2), base(3), R2(1), R2(2), R2(3), 'blue', 'LineWidth', 2)
        quiver3(base(1), base(2), base(3), R3(1), R3(2), R3(3), 'green', 'LineWidth', 2)
    end
end