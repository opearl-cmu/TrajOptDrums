function [YB_prescribed] = Prescribe_YB(t)
% Prescribe the kinematic trajectory
global t_YB_des;
global YB_des;
YB_prescribed = interp1(t_YB_des,YB_des,t,'spline','extrap');
end