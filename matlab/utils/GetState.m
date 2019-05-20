function state = GetState(homoTransform)
    if size(homoTransform,1)~=4
        homoTransform = reshape(homoTransform,4,4);
    end
    [r,p,w]=bot_quat_to_roll_pitch_yaw(bot_matrix_to_quat(homoTransform(1:3,1:3)));
    t=homoTransform(1:3,4);
    
    state = [t' r p w];