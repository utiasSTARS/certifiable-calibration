function t = ec_get_t_opt(R,Mtt,Mrt)
%ec_get_t_opt 

r = [vec(R); 1];
t = -inv(Mtt)*Mrt*r;

end

