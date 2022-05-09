function [F_D_prescribed] = Prescribe_F_D(t)
% Prescribe the kinematic trajectory
global t_F_D_des;
global F_D_des;
F_D_prescribed = interp1(t_F_D_des,F_D_des,t,'linear','extrap');
end