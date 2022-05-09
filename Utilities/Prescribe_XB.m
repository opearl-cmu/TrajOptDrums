function [XB_prescribed] = Prescribe_XB(t)
% Prescribe the kinematic trajectory
global t_XB_des;
global XB_des;
XB_prescribed = interp1(t_XB_des,XB_des,t,'spline','extrap');
end