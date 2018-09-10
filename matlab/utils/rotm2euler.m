function rpy = rotm2euler(R)
%rotm2euler 

[r,p,y]=bot_quat_to_roll_pitch_yaw(bot_matrix_to_quat(R));

rpy = [r p y];
end

