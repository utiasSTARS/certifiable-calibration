function R = euler2rotm(rpy)
%rotm2euler 
R = bot_quat_to_matrix(bot_roll_pitch_yaw_to_quat(rpy));
end

