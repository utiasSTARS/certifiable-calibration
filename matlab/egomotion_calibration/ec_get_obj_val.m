function val = ec_get_obj_val(R_calib,t_calib,R1,t1,R2,t2)
% Gets the value of the objective function as derived geometrically.

val = 0;

for i = 1:1:size(t1,2)
    temp = [R_calib*R2(:,:,i), R_calib*t2(:,i)+t_calib;zeros(1,3),1] - ...
        [R1(:,:,i)*R_calib, R1(:,:,i)*t_calib+t1(:,i);zeros(1,3),1];
    val = val + norm(temp,'fro');
end

end

