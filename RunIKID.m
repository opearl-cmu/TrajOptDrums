%% INVERSE KINEMATICS AND INVERSE DYNAMICS ON EXPERIMENTAL DATA
% Owen Douglas Pearl | CMU MBL | SPRING 2022
close all; clear all; clc

% Load in raw joint center data
load('SingleDownStroke_100bpm_JointCenters.mat')

% Put together experimental joint center table and smooth raw data
jointcenter_table = [smooth(elbow_x_exp, 10) smooth(elbow_y_exp, 10) smooth(wrist_x_exp, 10) smooth(wrist_y_exp, 10) smooth(finger_x_exp, 10) smooth(finger_y_exp, 10) smooth(bead_x_exp, 10) smooth(bead_y_exp, 10)] / 100;

% Set weights on kinematics
w = [1, 1, 1, 1]; 

% Setup optimization problem
optfunc = @(kins) ForwardKinematicsLoss(kins, time_exp, w, jointcenter_table);
kins_0 = zeros(4, length(time_exp));
options = optimoptions('fmincon', 'Display', 'off', 'MaxIter', 1000, 'MaxFunEvals',1e6, 'PlotFcns','optimplotfval');
[optimal_kinematics, optJ] = fmincon(optfunc, kins_0, [], [], [], [], [], [], [], options);

% COMPUTE KINEMATICS, DYNAMICS, AND COST & VISUALIZE
% Inputs for State Vizualization
t_viz_state = time_exp;
qs_viz_state = optimal_kinematics(1,:)';
qe_viz_state = optimal_kinematics(2,:)';
qw_viz_state = optimal_kinematics(3,:)';
qf_viz_state = optimal_kinematics(4,:)';
qsp_viz_state = gradient(qs_viz_state,t_viz_state);
qep_viz_state = gradient(qe_viz_state,t_viz_state);
qwp_viz_state = gradient(qw_viz_state,t_viz_state);
qfp_viz_state = gradient(qf_viz_state,t_viz_state);
qspp_viz_state = gradient(qsp_viz_state, t_viz_state);
qepp_viz_state = gradient(qep_viz_state, t_viz_state);
qwpp_viz_state = gradient(qwp_viz_state, t_viz_state);
qfpp_viz_state = gradient(qfp_viz_state, t_viz_state);
F_d_viz_state = zeros(size(t_viz_state));

% Use ID equations to calculate torque profiles from kinematics
SetParams
T_s_viz_state = 0.5.*G.*L_U.*M_U.*cos(qs_viz_state) + 0.5.*G.*M_L.*(2.*L_U.*cos(qs_viz_state)+L_L.*cos(qe_viz_state+qs_viz_state)) + 0.5.*G.*M_H.*(2.*L_U.*cos(qs_viz_state)+2.*L_L.*cos(qe_viz_state+qs_viz_state)+L_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state)) + 0.1666667.*F_d_viz_state.*(6.*L_L.*cos(qe_viz_state+qs_viz_state)+L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state)+3.*(L_S+2.*L_H).*cos(qe_viz_state+qs_viz_state+qw_viz_state)) + 0.1666667.*G.*M_S.*(6.*L_U.*cos(qs_viz_state)+6.*L_L.*cos(qe_viz_state+qs_viz_state)+6.*L_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state)+L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state)) + 0.5.*L_L.*L_U.*M_L.*sin(qe_viz_state).*(qsp_viz_state.^2-(qsp_viz_state+qep_viz_state).^2) + 0.5.*M_H.*(2.*L_L.*L_U.*sin(qe_viz_state).*qsp_viz_state.^2+L_H.*L_U.*sin(qe_viz_state+qw_viz_state).*qsp_viz_state.^2+L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state).^2-2.*L_L.*L_U.*sin(qe_viz_state).*(qsp_viz_state+qep_viz_state).^2-L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2-L_H.*L_U.*sin(qe_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2) + 0.1666667.*M_S.*(6.*L_L.*L_U.*sin(qe_viz_state).*qsp_viz_state.^2+6.*L_H.*L_U.*sin(qe_viz_state+qw_viz_state).*qsp_viz_state.^2+L_S.*L_U.*sin(qe_viz_state+qf_viz_state+qw_viz_state).*qsp_viz_state.^2+6.*L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state).^2+L_H.*L_S.*sin(qf_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2+L_L.*L_S.*sin(qf_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state).^2-6.*L_L.*L_U.*sin(qe_viz_state).*(qsp_viz_state+qep_viz_state).^2-6.*L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2-6.*L_H.*L_U.*sin(qe_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2-L_H.*L_S.*sin(qf_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state+qfp_viz_state).^2-L_L.*L_S.*sin(qf_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state+qfp_viz_state).^2-L_S.*L_U.*sin(qe_viz_state+qf_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state+qfp_viz_state).^2) + 0.02777778.*(36.*I_S+L_S.*M_S.*(L_S+6.*L_H.*cos(qf_viz_state)+6.*L_L.*cos(qf_viz_state+qw_viz_state)+6.*L_U.*cos(qe_viz_state+qf_viz_state+qw_viz_state))).*qwpp_viz_state + 0.02777778.*(36.*I_H+36.*I_S+9.*L_H.*M_H.*(L_H+2.*L_L.*cos(qw_viz_state)+2.*L_U.*cos(qe_viz_state+qw_viz_state))+M_S.*(L_S.^2+36.*L_H.^2+12.*L_H.*L_S.*cos(qf_viz_state)+36.*L_H.*L_L.*cos(qw_viz_state)+6.*L_L.*L_S.*cos(qf_viz_state+qw_viz_state)+36.*L_H.*L_U.*cos(qe_viz_state+qw_viz_state)+6.*L_S.*L_U.*cos(qe_viz_state+qf_viz_state+qw_viz_state))).*qwpp_viz_state + 0.02777778.*(36.*I_H+36.*I_L+36.*I_S+9.*L_L.*M_L.*(L_L+2.*L_U.*cos(qe_viz_state))+9.*M_H.*(L_H.^2+4.*L_L.^2+4.*L_H.*L_L.*cos(qw_viz_state)+4.*L_L.*L_U.*cos(qe_viz_state)+2.*L_H.*L_U.*cos(qe_viz_state+qw_viz_state))+M_S.*(L_S.^2+36.*L_H.^2+36.*L_L.^2+12.*L_H.*L_S.*cos(qf_viz_state)+36.*L_L.*L_U.*cos(qe_viz_state)+72.*L_H.*L_L.*cos(qw_viz_state)+12.*L_L.*L_S.*cos(qf_viz_state+qw_viz_state)+36.*L_H.*L_U.*cos(qe_viz_state+qw_viz_state)+6.*L_S.*L_U.*cos(qe_viz_state+qf_viz_state+qw_viz_state))).*qepp_viz_state + 0.02777778.*(36.*I_H+36.*I_L+36.*I_S+36.*I_U+9.*M_U.*L_U.^2+9.*M_L.*(L_L.^2+4.*L_U.^2+4.*L_L.*L_U.*cos(qe_viz_state))+9.*M_H.*(L_H.^2+4.*L_L.^2+4.*L_U.^2+4.*L_H.*L_L.*cos(qw_viz_state)+8.*L_L.*L_U.*cos(qe_viz_state)+4.*L_H.*L_U.*cos(qe_viz_state+qw_viz_state))+M_S.*(L_S.^2+36.*L_H.^2+36.*L_L.^2+36.*L_U.^2+12.*L_H.*L_S.*cos(qf_viz_state)+72.*L_H.*L_L.*cos(qw_viz_state)+72.*L_L.*L_U.*cos(qe_viz_state)+12.*L_L.*L_S.*cos(qf_viz_state+qw_viz_state)+72.*L_H.*L_U.*cos(qe_viz_state+qw_viz_state)+12.*L_S.*L_U.*cos(qe_viz_state+qf_viz_state+qw_viz_state))).*qspp_viz_state - 0.5.*G.*L_L.*M_L.*cos(qe_viz_state+qs_viz_state) - 0.5.*G.*M_H.*(2.*L_L.*cos(qe_viz_state+qs_viz_state)+L_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state)) - 0.1666667.*G.*M_S.*(6.*L_L.*cos(qe_viz_state+qs_viz_state)+6.*L_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state)+L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state)) - 0.1666667.*F_d_viz_state.*(6.*L_U.*cos(qs_viz_state)+6.*L_L.*cos(qe_viz_state+qs_viz_state)+L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state)+3.*(L_S+2.*L_H).*cos(qe_viz_state+qs_viz_state+qw_viz_state)) - 0.5.*L_L.*L_U.*M_L.*sin(qe_viz_state).*qsp_viz_state.^2 - 0.5.*M_H.*(2.*L_L.*L_U.*sin(qe_viz_state).*qsp_viz_state.^2+L_H.*L_U.*sin(qe_viz_state+qw_viz_state).*qsp_viz_state.^2+L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state).^2-L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2) - 0.1666667.*M_S.*(6.*L_L.*L_U.*sin(qe_viz_state).*qsp_viz_state.^2+6.*L_H.*L_U.*sin(qe_viz_state+qw_viz_state).*qsp_viz_state.^2+L_S.*L_U.*sin(qe_viz_state+qf_viz_state+qw_viz_state).*qsp_viz_state.^2+6.*L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state).^2+L_H.*L_S.*sin(qf_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2+L_L.*L_S.*sin(qf_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state).^2-6.*L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2-L_H.*L_S.*sin(qf_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state+qfp_viz_state).^2-L_L.*L_S.*sin(qf_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state+qfp_viz_state).^2) - 0.02777778.*(36.*I_S+L_S.*M_S.*(L_S+6.*L_H.*cos(qf_viz_state)+6.*L_L.*cos(qf_viz_state+qw_viz_state))).*qwpp_viz_state - 0.02777778.*(36.*I_H+36.*I_S+9.*L_H.*M_H.*(L_H+2.*L_L.*cos(qw_viz_state))+M_S.*(L_S.^2+36.*L_H.^2+12.*L_H.*L_S.*cos(qf_viz_state)+36.*L_H.*L_L.*cos(qw_viz_state)+6.*L_L.*L_S.*cos(qf_viz_state+qw_viz_state))).*qwpp_viz_state - 0.02777778.*(36.*I_H+36.*I_L+36.*I_S+9.*M_L.*L_L.^2+9.*M_H.*(L_H.^2+4.*L_L.^2+4.*L_H.*L_L.*cos(qw_viz_state))+M_S.*(L_S.^2+36.*L_H.^2+36.*L_L.^2+12.*L_H.*L_S.*cos(qf_viz_state)+72.*L_H.*L_L.*cos(qw_viz_state)+12.*L_L.*L_S.*cos(qf_viz_state+qw_viz_state))).*qepp_viz_state - 0.02777778.*(36.*I_H+36.*I_L+36.*I_S+9.*L_L.*M_L.*(L_L+2.*L_U.*cos(qe_viz_state))+9.*M_H.*(L_H.^2+4.*L_L.^2+4.*L_H.*L_L.*cos(qw_viz_state)+4.*L_L.*L_U.*cos(qe_viz_state)+2.*L_H.*L_U.*cos(qe_viz_state+qw_viz_state))+M_S.*(L_S.^2+36.*L_H.^2+36.*L_L.^2+12.*L_H.*L_S.*cos(qf_viz_state)+36.*L_L.*L_U.*cos(qe_viz_state)+72.*L_H.*L_L.*cos(qw_viz_state)+12.*L_L.*L_S.*cos(qf_viz_state+qw_viz_state)+36.*L_H.*L_U.*cos(qe_viz_state+qw_viz_state)+6.*L_S.*L_U.*cos(qe_viz_state+qf_viz_state+qw_viz_state))).*qspp_viz_state;
T_e_viz_state = 0.5.*G.*L_L.*M_L.*cos(qe_viz_state+qs_viz_state) + 0.5.*G.*M_H.*(2.*L_L.*cos(qe_viz_state+qs_viz_state)+L_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state)) + 0.1666667.*F_d_viz_state.*(L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state)+3.*(L_S+2.*L_H).*cos(qe_viz_state+qs_viz_state+qw_viz_state)) + 0.1666667.*G.*M_S.*(6.*L_L.*cos(qe_viz_state+qs_viz_state)+6.*L_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state)+L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state)) + 0.5.*L_L.*L_U.*M_L.*sin(qe_viz_state).*qsp_viz_state.^2 + 0.5.*M_H.*(2.*L_L.*L_U.*sin(qe_viz_state).*qsp_viz_state.^2+L_H.*L_U.*sin(qe_viz_state+qw_viz_state).*qsp_viz_state.^2+L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state).^2-L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2) + 0.1666667.*M_S.*(6.*L_L.*L_U.*sin(qe_viz_state).*qsp_viz_state.^2+6.*L_H.*L_U.*sin(qe_viz_state+qw_viz_state).*qsp_viz_state.^2+L_S.*L_U.*sin(qe_viz_state+qf_viz_state+qw_viz_state).*qsp_viz_state.^2+6.*L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state).^2+L_H.*L_S.*sin(qf_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2+L_L.*L_S.*sin(qf_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state).^2-6.*L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2-L_H.*L_S.*sin(qf_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state+qfp_viz_state).^2-L_L.*L_S.*sin(qf_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state+qfp_viz_state).^2) + 0.02777778.*(36.*I_S+L_S.*M_S.*(L_S+6.*L_H.*cos(qf_viz_state)+6.*L_L.*cos(qf_viz_state+qw_viz_state))).*qwpp_viz_state + 0.02777778.*(36.*I_H+36.*I_S+9.*L_H.*M_H.*(L_H+2.*L_L.*cos(qw_viz_state))+M_S.*(L_S.^2+36.*L_H.^2+12.*L_H.*L_S.*cos(qf_viz_state)+36.*L_H.*L_L.*cos(qw_viz_state)+6.*L_L.*L_S.*cos(qf_viz_state+qw_viz_state))).*qwpp_viz_state + 0.02777778.*(36.*I_H+36.*I_L+36.*I_S+9.*M_L.*L_L.^2+9.*M_H.*(L_H.^2+4.*L_L.^2+4.*L_H.*L_L.*cos(qw_viz_state))+M_S.*(L_S.^2+36.*L_H.^2+36.*L_L.^2+12.*L_H.*L_S.*cos(qf_viz_state)+72.*L_H.*L_L.*cos(qw_viz_state)+12.*L_L.*L_S.*cos(qf_viz_state+qw_viz_state))).*qepp_viz_state + 0.02777778.*(36.*I_H+36.*I_L+36.*I_S+9.*L_L.*M_L.*(L_L+2.*L_U.*cos(qe_viz_state))+9.*M_H.*(L_H.^2+4.*L_L.^2+4.*L_H.*L_L.*cos(qw_viz_state)+4.*L_L.*L_U.*cos(qe_viz_state)+2.*L_H.*L_U.*cos(qe_viz_state+qw_viz_state))+M_S.*(L_S.^2+36.*L_H.^2+36.*L_L.^2+12.*L_H.*L_S.*cos(qf_viz_state)+36.*L_L.*L_U.*cos(qe_viz_state)+72.*L_H.*L_L.*cos(qw_viz_state)+12.*L_L.*L_S.*cos(qf_viz_state+qw_viz_state)+36.*L_H.*L_U.*cos(qe_viz_state+qw_viz_state)+6.*L_S.*L_U.*cos(qe_viz_state+qf_viz_state+qw_viz_state))).*qspp_viz_state - 0.5.*G.*L_H.*M_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state) - 0.1666667.*G.*M_S.*(6.*L_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state)+L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state)) - 0.1666667.*F_d_viz_state.*(6.*L_L.*cos(qe_viz_state+qs_viz_state)+L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state)+3.*(L_S+2.*L_H).*cos(qe_viz_state+qs_viz_state+qw_viz_state)) - 0.5.*L_H.*M_H.*(L_U.*sin(qe_viz_state+qw_viz_state).*qsp_viz_state.^2+L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state).^2) - 0.1666667.*M_S.*(6.*L_H.*L_U.*sin(qe_viz_state+qw_viz_state).*qsp_viz_state.^2+L_S.*L_U.*sin(qe_viz_state+qf_viz_state+qw_viz_state).*qsp_viz_state.^2+6.*L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state).^2+L_H.*L_S.*sin(qf_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2+L_L.*L_S.*sin(qf_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state).^2-L_H.*L_S.*sin(qf_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state+qfp_viz_state).^2) - 0.02777778.*(36.*I_S+L_S.*M_S.*(L_S+6.*L_H.*cos(qf_viz_state))).*qwpp_viz_state - 0.02777778.*(36.*I_H+36.*I_S+9.*M_H.*L_H.^2+M_S.*(L_S.^2+36.*L_H.^2+12.*L_H.*L_S.*cos(qf_viz_state))).*qwpp_viz_state - 0.02777778.*(36.*I_H+36.*I_S+9.*L_H.*M_H.*(L_H+2.*L_L.*cos(qw_viz_state))+M_S.*(L_S.^2+36.*L_H.^2+12.*L_H.*L_S.*cos(qf_viz_state)+36.*L_H.*L_L.*cos(qw_viz_state)+6.*L_L.*L_S.*cos(qf_viz_state+qw_viz_state))).*qepp_viz_state - 0.02777778.*(36.*I_H+36.*I_S+9.*L_H.*M_H.*(L_H+2.*L_L.*cos(qw_viz_state)+2.*L_U.*cos(qe_viz_state+qw_viz_state))+M_S.*(L_S.^2+36.*L_H.^2+12.*L_H.*L_S.*cos(qf_viz_state)+36.*L_H.*L_L.*cos(qw_viz_state)+6.*L_L.*L_S.*cos(qf_viz_state+qw_viz_state)+36.*L_H.*L_U.*cos(qe_viz_state+qw_viz_state)+6.*L_S.*L_U.*cos(qe_viz_state+qf_viz_state+qw_viz_state))).*qspp_viz_state;
T_w_viz_state = 0.5.*G.*L_H.*M_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state) + 0.1666667.*F_d_viz_state.*L_S.*(cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state)+3.*cos(qe_viz_state+qs_viz_state+qw_viz_state)) + 0.1666667.*G.*M_S.*(6.*L_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state)+L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state)) + 0.5.*L_H.*M_H.*(L_U.*sin(qe_viz_state+qw_viz_state).*qsp_viz_state.^2+L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state).^2) + 0.1666667.*M_S.*(6.*L_H.*L_U.*sin(qe_viz_state+qw_viz_state).*qsp_viz_state.^2+L_S.*L_U.*sin(qe_viz_state+qf_viz_state+qw_viz_state).*qsp_viz_state.^2+6.*L_H.*L_L.*sin(qw_viz_state).*(qsp_viz_state+qep_viz_state).^2+L_H.*L_S.*sin(qf_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2+L_L.*L_S.*sin(qf_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state).^2-L_H.*L_S.*sin(qf_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state+qfp_viz_state).^2) + 0.02777778.*(36.*I_S+L_S.*M_S.*(L_S+6.*L_H.*cos(qf_viz_state))).*qwpp_viz_state + 0.02777778.*(36.*I_H+36.*I_S+9.*M_H.*L_H.^2+M_S.*(L_S.^2+36.*L_H.^2+12.*L_H.*L_S.*cos(qf_viz_state))).*qwpp_viz_state + 0.02777778.*(36.*I_H+36.*I_S+9.*L_H.*M_H.*(L_H+2.*L_L.*cos(qw_viz_state))+M_S.*(L_S.^2+36.*L_H.^2+12.*L_H.*L_S.*cos(qf_viz_state)+36.*L_H.*L_L.*cos(qw_viz_state)+6.*L_L.*L_S.*cos(qf_viz_state+qw_viz_state))).*qepp_viz_state + 0.02777778.*(36.*I_H+36.*I_S+9.*L_H.*M_H.*(L_H+2.*L_L.*cos(qw_viz_state)+2.*L_U.*cos(qe_viz_state+qw_viz_state))+M_S.*(L_S.^2+36.*L_H.^2+12.*L_H.*L_S.*cos(qf_viz_state)+36.*L_H.*L_L.*cos(qw_viz_state)+6.*L_L.*L_S.*cos(qf_viz_state+qw_viz_state)+36.*L_H.*L_U.*cos(qe_viz_state+qw_viz_state)+6.*L_S.*L_U.*cos(qe_viz_state+qf_viz_state+qw_viz_state))).*qspp_viz_state - 0.1666667.*G.*L_S.*M_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state) - 0.1666667.*F_d_viz_state.*(L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state)+3.*(L_S+2.*L_H).*cos(qe_viz_state+qs_viz_state+qw_viz_state)) - 0.1666667.*L_S.*M_S.*(L_U.*sin(qe_viz_state+qf_viz_state+qw_viz_state).*qsp_viz_state.^2+L_H.*sin(qf_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2+L_L.*sin(qf_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state).^2) - 0.02777778.*(36.*I_S+M_S.*L_S.^2).*qwpp_viz_state - 0.02777778.*(36.*I_S+L_S.*M_S.*(L_S+6.*L_H.*cos(qf_viz_state))).*qwpp_viz_state - 0.02777778.*(36.*I_S+L_S.*M_S.*(L_S+6.*L_H.*cos(qf_viz_state)+6.*L_L.*cos(qf_viz_state+qw_viz_state))).*qepp_viz_state - 0.02777778.*(36.*I_S+L_S.*M_S.*(L_S+6.*L_H.*cos(qf_viz_state)+6.*L_L.*cos(qf_viz_state+qw_viz_state)+6.*L_U.*cos(qe_viz_state+qf_viz_state+qw_viz_state))).*qspp_viz_state;
T_f_viz_state = 0.1666667.*G.*L_S.*M_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state) + 0.1666667.*L_S.*M_S.*(L_U.*sin(qe_viz_state+qf_viz_state+qw_viz_state).*qsp_viz_state.^2+L_H.*sin(qf_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state).^2+L_L.*sin(qf_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state).^2) + 0.02777778.*(36.*I_S+M_S.*L_S.^2).*qwpp_viz_state + 0.02777778.*(36.*I_S+L_S.*M_S.*(L_S+6.*L_H.*cos(qf_viz_state))).*qwpp_viz_state + 0.02777778.*(36.*I_S+L_S.*M_S.*(L_S+6.*L_H.*cos(qf_viz_state)+6.*L_L.*cos(qf_viz_state+qw_viz_state))).*qepp_viz_state + 0.02777778.*(36.*I_S+L_S.*M_S.*(L_S+6.*L_H.*cos(qf_viz_state)+6.*L_L.*cos(qf_viz_state+qw_viz_state)+6.*L_U.*cos(qe_viz_state+qf_viz_state+qw_viz_state))).*qspp_viz_state - 0.1666667.*F_d_viz_state.*L_S.*(cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state)+3.*cos(qe_viz_state+qs_viz_state+qw_viz_state));

% Calculate Drum Bead XY-Pose
xb_viz_state = L_U.*cos(qs_viz_state) + L_L.*cos(qe_viz_state+qs_viz_state) + L_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state) + 0.6666667.*L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state);
yb_viz_state = L_U.*sin(qs_viz_state) + L_L.*sin(qe_viz_state+qs_viz_state) + L_H.*sin(qe_viz_state+qs_viz_state+qw_viz_state) + 0.6666667.*L_S.*sin(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state);

% Calculate butt of drumstick pose
xbutt_viz_state = L_U.*cos(qs_viz_state) + L_L.*cos(qe_viz_state+qs_viz_state) + L_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state) - 0.3333333.*L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state);
ybutt_viz_state = L_U.*sin(qs_viz_state) + L_L.*sin(qe_viz_state+qs_viz_state) + L_H.*sin(qe_viz_state+qs_viz_state+qw_viz_state) - 0.3333333.*L_S.*sin(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state);

% Calculate all joint center positions for vizualization
xe_viz_state = L_U.*cos(qs_viz_state);
ye_viz_state = L_U.*sin(qs_viz_state);
xw_viz_state = L_U.*cos(qs_viz_state) + L_L.*cos(qe_viz_state + qs_viz_state);
yw_viz_state = L_U.*sin(qs_viz_state) + L_L.*sin(qe_viz_state + qs_viz_state);
xf_viz_state = L_U.*cos(qs_viz_state) + L_L.*cos(qe_viz_state + qs_viz_state) + L_H.*cos(qe_viz_state + qs_viz_state + qw_viz_state);
yf_viz_state = L_U.*sin(qs_viz_state) + L_L.*sin(qe_viz_state + qs_viz_state) + L_H.*sin(qe_viz_state + qs_viz_state + qw_viz_state);

% Calculate total (integral) of effort cost (sum of square torques)
Jeffort = smooth(T_s_viz_state, 10).^2 + smooth(T_e_viz_state, 10).^2 + smooth(T_w_viz_state, 10).^2 + smooth(T_f_viz_state, 10).^2;

% Visualize Stroke Trajectories
lim_y = [-.6 0.1]; lim_x = [-.2 .8]; frames = [130];
VisualizeDrumStroke

disp('IK & ID Completed')

%% DEFINE FORWARD KINEMATICS LOSS FUNCTION
function [J] = ForwardKinematicsLoss(kins, t, w, jointcenter_table)
% Load Model Parameters
SetParams

% Unpack kins and differentiate for velocities/accelerations
qs = kins(1,:);
qe = kins(2,:);
qw = kins(3,:);
qf = kins(4,:);

% Calculate all joint center positions with kins model from given kins
xe = L_U.*cos(qs);
ye = L_U.*sin(qs);
xw = L_U.*cos(qs) + L_L.*cos(qe + qs);
yw = L_U.*sin(qs) + L_L.*sin(qe + qs);
xf = L_U.*cos(qs) + L_L.*cos(qe + qs) + L_H.*cos(qe + qs + qw);
yf = L_U.*sin(qs) + L_L.*sin(qe + qs) + L_H.*sin(qe + qs + qw);
xb = L_U.*cos(qs) + L_L.*cos(qe+qs) + L_H.*cos(qe+qs+qw) + 0.6666667.*L_S.*cos(qe+qf+qs+qw);
yb = L_U.*sin(qs) + L_L.*sin(qe+qs) + L_H.*sin(qe+qs+qw) + 0.6666667.*L_S.*sin(qe+qf+qs+qw);

% Unpack jointcenter_table
elbow_x_exp = jointcenter_table(:,1);
elbow_y_exp = jointcenter_table(:,2);
wrist_x_exp = jointcenter_table(:,3);
wrist_y_exp = jointcenter_table(:,4);
finger_x_exp = jointcenter_table(:,5);
finger_y_exp = jointcenter_table(:,6);
bead_x_exp = jointcenter_table(:,7);
bead_y_exp = jointcenter_table(:,8);

% Calculate loss
Cost = w(1).*( (elbow_x_exp - xe').^2 + (elbow_y_exp - ye').^2 ) + w(2).*( (wrist_x_exp - xw').^2 + (wrist_y_exp - yw').^2 ) + w(3).*( (finger_x_exp - xf').^2 + (finger_y_exp - yf').^2 ) + w(4).*( (bead_x_exp - xb').^2 + (bead_y_exp - yb').^2 );
J = trapz(t, Cost);
end