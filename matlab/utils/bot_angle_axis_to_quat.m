function q = bot_angle_axis_to_quat (theta, axis)
    x = axis(1+0);
    y = axis(1+1); 
    z = axis(1+2);
    norm = sqrt (x*x + y*y + z*z);
    if (0 == norm)
        q(1+0) = 1;
        q(1+1) = 0;
        q(1+2) = 0;
        q(1+3) = 0;
        return;
    end

    t = sin(theta/2) / norm;

    q(1+0) = cos(theta / 2);
    q(1+1) = x * t;
    q(1+2) = y * t;
    q(1+3) = z * t;