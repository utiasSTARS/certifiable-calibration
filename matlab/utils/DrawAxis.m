function DrawAxis(T, scale, cx, cy, cz)    

    origin = T * [0; 0; 0; 1];
    xaxis = T * [scale; 0; 0; 1];
    yaxis = T * [0; scale; 0; 1];
    zaxis = T * [0; 0; scale; 1];

    a = [origin xaxis];
    line(a(1,:), a(2,:), a(3,:), 'Color', cx, 'LineWidth', 3);
    a = [origin yaxis];
    line(a(1,:), a(2,:), a(3,:), 'Color', cy, 'LineWidth', 3);
    a = [origin zaxis];
    line(a(1,:), a(2,:), a(3,:), 'Color', cz, 'LineWidth', 3);
