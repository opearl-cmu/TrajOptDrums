function [Cost, Cost_Grad] = DrummingArmCostFunc(t, z, u)
% Unpack Controls
T_S  = u(1,:);
T_E = u(2,:);
T_W = u(3,:);
T_F = u(4,:);

% Calculate the Cost Given the Current State (min efforrt)
Cost = (T_S).^2 + (T_E).^2 + (T_W).^2 + (T_F).^2;

% Calculate the Cost Gradient
CalculateCostGradients;
Cost_Grad = [DJDT; DJDQS; DJDQE; DJDQW; DJDQF; DJDQSP; DJDQEP; DJDQWP; DJDQFP; DJDT_S; DJDT_E; DJDT_W; DJDT_F; DJDF_D];
end