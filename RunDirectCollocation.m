%% DRUMMING ARM BIOMECHANICAL MODEL TRAJECTORY OPTIMIZATION
% Owen Douglas Pearl | CMU MBL | SPRING 2022
close all; clear all; clc

% Load experimental trajectory for drumstick bead for tracking
load('ExperimentalData\ExpIKBeadTrajectory.mat');
global t_YB_des; t_YB_des = t_grid;
global YB_des; YB_des = smooth(yb_grid);
global t_F_D_des; t_F_D_des = t_grid;
global F_D_des; F_D_des = zeros(size(t_grid));
global t_XB_des; t_XB_des = t_grid;
global XB_des; XB_des = smooth(xb_grid);
t0 = 0;  tF = t_grid(end); ttol = 0.0;
qs0 = qs_grid(1);
qe0 = qe_grid(1);
qw0 = qw_grid(1);
qtol = 0;
qf0 = qf_grid(1);
qsp0 = qsp_grid(1);
qep0 = qep_grid(1);
qwp0 = qwp_grid(1);
qfp0 = qfp_grid(1);
qptol = 0;
clearvars -except t_YB_des YB_des XB_des t_F_D_des F_D_des t0 tF ttol qs0 qe0 qw0 qtol qf0 qsp0 qep0 qwp0 qfp0 qptol

% profile OFF

% Run Direct Collocation (using OptimTraj library)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%             Set up optimization params and lookup tables                %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Load Model Parameters
SetParams

% Track XB as well as YB?
trackXB = true;

% Create output plots/vids at the end of simulation
createPlots = 1;

% Save output of optimization as .mat
doSave = 0;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set up function handles                           %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
problem.func.dynamics = @(t,z,u)( DrummingArmDynamics(t, z, u) );

problem.func.pathObj = @(t,z,u)( DrummingArmCostFunc(t, z, u) );

problem.func.pathCst = @(t,z,u)( DrummingArmConstraints(t, z, u, trackXB) );

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%               Set up bounds on time, state, and control                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
problem.bounds.initialTime.low = t0;
problem.bounds.initialTime.upp = t0 + ttol;
problem.bounds.finalTime.low = tF - ttol;
problem.bounds.finalTime.upp = tF + ttol;

problem.bounds.state.low = [qs0 - pi/6; qe0 - pi/6; 0.0017; qf0; -3; -8; -17; -23];
problem.bounds.state.upp = [qs0 + pi/6; qe0 + pi/3; qw0 + pi/4; qf0 + pi/6; 3; 8; 17; 23];

problem.bounds.initialState.low = [qs0 - qtol; qe0 - qtol; qw0 - qtol; qf0 - qtol; qsp0 - qptol; qep0 - qptol; qwp0 - qptol; qfp0 - qptol];
problem.bounds.initialState.upp = [qs0 + qtol; qe0 + qtol; qw0  + qtol; qf0 + qtol; qsp0 + qptol; qep0 + qptol; qwp0 + qptol; qfp0 + qptol];

T_S_low = -10; T_S_upp = 10;
T_E_low = -10; T_E_upp = 10;
T_W_low = -10; T_W_upp = 10;
T_F_low = -5; T_F_upp = 0;
F_D_low = 0; F_D_upp = inf; 
problem.bounds.control.low = [T_S_low; T_E_low; T_W_low; T_F_low; F_D_low];
problem.bounds.control.upp = [T_S_upp; T_E_upp; T_W_upp; T_F_upp; F_D_upp];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%              Create an initial guess for the trajectory                 %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Most basic possible initial guess 
% t_guess = [t0, tF];
% z_guess = [qs0, qs0; qe0, qe0; qw0, qw0; qf0, qf0; qsp0, qsp0; qep0, qep0; qwp0, qwp0; qfp0, qfp0];
% u_guess = [0, 0; 0, 0; 0, 0; 0, 0; 0, 0];

% Warmstart with exp solution from TVLQR
load('ExperimentalData\TVLQR_Result.mat', 'soln');
t_guess = soln.grid.time';
z_guess = soln.grid.state;
u_guess = soln.grid.control;
clearvars 'soln'

problem.guess.time = t_guess;
problem.guess.state = z_guess;
problem.guess.control = u_guess;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Options:                                      %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% method = 'trapGrad';
method = 'hermiteSimpsonGrad'; segN1 = 15; segN2 = 100;
% method = 'rungeKuttaGrad';

%%% Decide whether or not to use IPOPT (if fmincon; use less segments)
problem.options(1).useIPOPT = false;
problem.options(2).useIPOPT = true; problem.options(2).trackXB = trackXB;

%%% Method-independent options:
problem.options(1).nlpOpt = optimset(...
    'Display','iter',...   % {'iter','final','off'}
    'TolFun',1e-3,...
    'MaxFunEvals',1e4,...
    'MaxIter', 500,...
     'TolCon', 3); 
%     'PlotFcns','optimplotfval');   %options for fmincon
problem.options(1).verbose = 3; % How much to print out?

problem.options(2).nlpOpt = optimset(...
    'Display','iter',...   % {'iter','final','off'}
    'TolFun',1e-4,...
    'MaxFunEvals',5e6,...   %options for fmincon
    'MaxIter', 3000); % 
%     'PlotFcns','optimplotfval');
problem.options(2).verbose = 3; % How much to print out

switch method
        
    case 'trapGrad'  %trapezoid with analytic gradients
        
        problem.options(1).method = 'trapezoid'; % Select the transcription method
        problem.options(1).trapezoid.nGrid = 10;  %method-specific options
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        
        problem.options(2).method = 'trapezoid'; % Select the transcription method
        problem.options(2).trapezoid.nGrid = 45;  %method-specific options
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';
        
    case 'hermiteSimpsonGrad'  %hermite simpson with analytic gradients
        
        problem.options(1).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(1).hermiteSimpson.nSegment = segN1;  %method-specific options
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        
        problem.options(2).method = 'hermiteSimpson'; % Select the transcription method
        problem.options(2).hermiteSimpson.nSegment = segN2;  %method-specific options
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';
        problem.options(2).nlpOpt.DerivativeCheck = 'off';
    
    case 'rungeKuttaGrad'
      
        problem.options(1).method = 'rungeKutta'; % Select the transcription method
        problem.options(1).defaultAccuracy = 'low';
        problem.options(1).nlpOpt.GradConstr = 'on';
        problem.options(1).nlpOpt.GradObj = 'on';
        problem.options(1).nlpOpt.DerivativeCheck = 'off';
        
        problem.options(2).method = 'rungeKutta'; % Select the transcription method
        problem.options(2).defaultAccuracy = 'medium';
        problem.options(2).nlpOpt.GradConstr = 'on';
        problem.options(2).nlpOpt.GradObj = 'on';
        
    otherwise
        error('Invalid method!');
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Solve!                                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%%%%% THE KEY LINE:
tstart = tic;
soln = optimTraj(problem);
RunTime = toc(tstart);

if doSave == 1
    save('OptimizationResults.mat')
end

% profile VIEWER

%% Output, Compare, and Visualize Tracking Performance
% Transcription Grid points:
t_grid = soln(end).grid.time;
qs_grid = soln(end).grid.state(1,:);
qe_grid = soln(end).grid.state(2,:);
qw_grid = soln(end).grid.state(3,:);
qf_grid = soln(end).grid.state(4,:);
qsp_grid = soln(end).grid.state(5,:);
qep_grid = soln(end).grid.state(6,:);
qwp_grid = soln(end).grid.state(7,:);
qfp_grid = soln(end).grid.state(8,:);
T_s_grid = soln(end).grid.control(1,:);
T_e_grid = soln(end).grid.control(2,:);
T_w_grid = soln(end).grid.control(3,:);
T_f_grid = soln(end).grid.control(4,:);
F_d_grid = soln(end).grid.control(5,:);

% Inputs for State Vizualization
t_viz_state = t0:.01:tF;
qs_viz_state = smooth(interp1(t_grid, qs_grid, t_viz_state, 'spline', 'extrap'), 10);
qe_viz_state = smooth(interp1(t_grid, qe_grid, t_viz_state, 'spline', 'extrap'), 10);
qw_viz_state = smooth(interp1(t_grid, qw_grid, t_viz_state, 'spline', 'extrap'), 10);
qf_viz_state = smooth(interp1(t_grid, qf_grid, t_viz_state, 'spline', 'extrap'), 10);
qsp_viz_state = smooth(interp1(t_grid, qsp_grid, t_viz_state, 'spline', 'extrap'), 10);
qep_viz_state = smooth(interp1(t_grid, qep_grid, t_viz_state, 'spline', 'extrap'), 10);
qwp_viz_state = smooth(interp1(t_grid, qwp_grid, t_viz_state, 'spline', 'extrap'), 10);
qfp_viz_state = smooth(interp1(t_grid, qfp_grid, t_viz_state, 'spline', 'extrap'), 10);
T_s_viz_state = smooth(interp1(t_grid, T_s_grid, t_viz_state, 'spline', 'extrap'), 10);
T_e_viz_state = smooth(interp1(t_grid, T_e_grid, t_viz_state, 'spline', 'extrap'), 10);
T_w_viz_state = smooth(interp1(t_grid, T_w_grid, t_viz_state, 'spline', 'extrap'), 10);
T_f_viz_state = smooth(interp1(t_grid, T_f_grid, t_viz_state, 'spline', 'extrap'), 10);
F_d_viz_state = smooth(interp1(t_grid, F_d_grid, t_viz_state, 'spline', 'extrap'), 10);

% Calculate Drum Bead XY-Pose
xb_viz_state = L_U.*cos(qs_viz_state) + L_L.*cos(qe_viz_state+qs_viz_state) + L_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state) + 0.6666667.*L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state);
yb_viz_state = L_U.*sin(qs_viz_state) + L_L.*sin(qe_viz_state+qs_viz_state) + L_H.*sin(qe_viz_state+qs_viz_state+qw_viz_state) + 0.6666667.*L_S.*sin(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state);
ybp_viz_state = L_U.*cos(qs_viz_state).*qsp_viz_state + L_L.*cos(qe_viz_state+qs_viz_state).*(qsp_viz_state+qep_viz_state) + L_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state) + 0.6666667.*L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state).*(qsp_viz_state+qep_viz_state+qwp_viz_state+qfp_viz_state);
% ybpp_viz_state = -l_u.*(sin(qs).*qsp.^2-cos(qs).*qspp) - l_l.*(sin(qe+qs).*(qsp+qep).^2-cos(qe+qs).*(qepp+qspp)) - l_h.*(sin(qe+qs+qw).*(qsp+qep+qwp).^2-cos(qe+qs+qw).*(qepp+qspp+qwpp)) - 0.6666667.*l_s.*(sin(qe+qf+qs+qw).*(qsp+qep+qwp+qfp).^2-cos(qe+qf+qs+qw).*(qepp+qfpp+qspp+qwpp));

% Calculate prescribed drum bead trajectory
yb_viz_des = Prescribe_YB(t_viz_state);
% xb_viz_des = Prescribe_XB(t_viz_state);

% Calculate all joint center positions for simulation
xe_viz_state = L_U.*cos(qs_viz_state);
ye_viz_state = L_U.*sin(qs_viz_state);
xw_viz_state = L_U.*cos(qs_viz_state) + L_L.*cos(qe_viz_state + qs_viz_state);
yw_viz_state = L_U.*sin(qs_viz_state) + L_L.*sin(qe_viz_state + qs_viz_state);
xf_viz_state = L_U.*cos(qs_viz_state) + L_L.*cos(qe_viz_state + qs_viz_state) + L_H.*cos(qe_viz_state + qs_viz_state + qw_viz_state);
yf_viz_state = L_U.*sin(qs_viz_state) + L_L.*sin(qe_viz_state + qs_viz_state) + L_H.*sin(qe_viz_state + qs_viz_state + qw_viz_state);

% Calculate butt of drumstick pose
xbutt_viz_state = L_U.*cos(qs_viz_state) + L_L.*cos(qe_viz_state+qs_viz_state) + L_H.*cos(qe_viz_state+qs_viz_state+qw_viz_state) - 0.3333333.*L_S.*cos(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state);
ybutt_viz_state = L_U.*sin(qs_viz_state) + L_L.*sin(qe_viz_state+qs_viz_state) + L_H.*sin(qe_viz_state+qs_viz_state+qw_viz_state) - 0.3333333.*L_S.*sin(qe_viz_state+qf_viz_state+qs_viz_state+qw_viz_state);

% Calculate total (integral) of effort cost (sum of square torques)
Jeffort = smooth(T_s_viz_state, 10).^2 + smooth(T_e_viz_state, 10).^2 + smooth(T_f_viz_state, 10).^2 + smooth(T_w_viz_state, 10).^2;
Jrms = rms(Jeffort)

if(createPlots)
    % Plot drum bead trajectory
    figure 
    plot(smooth(xb_viz_state), smooth(yb_viz_state), 'k');
    title('Drum Stick Bead XY-Position')
    ylabel('Bead Y-Pose (y)')
    xlabel('Bead X-Pose (m)')
    ylim([-.5, .1])
    xlim([.4, .7])
    set(gca, 'Box', 'off')
    hold off
    
    % Plot Y-Trajectory against prescribed
    figure 
    plot(t_viz_state, smooth(yb_viz_state), 'k--');
    hold on
    plot(t_viz_state, smooth(yb_viz_des), 'k');
    ylabel('Bead Y-Pose Tracking (m)')
    xlabel('Time (s)')
    legend('State', 'Prescribed')
    set(gca, 'Box', 'off')
    legend boxoff
    hold off
    
    % Plot Kinematics
    figure
    RadToDeg = (180)/(pi);
    plot(t_viz_state, smooth(qs_viz_state.*RadToDeg, 10)+90, 'r')
    hold on
    plot(t_viz_state, smooth(qe_viz_state.*RadToDeg, 10)-90, 'g')
    plot(t_viz_state, smooth(qw_viz_state.*RadToDeg, 10), 'b')
    plot(t_viz_state, smooth(qf_viz_state.*RadToDeg, 10), 'k')
    legend('QS', 'QE', 'QW', 'QF')
    ylim([-50 60])
    ylabel('Joint Angles (deg)')
    xlabel('Time (s)')
    set(gca, 'Box', 'off')
    legend boxoff
    hold off
    
     % Plot Torques
    figure
    plot(t_viz_state, smooth(T_s_viz_state, 10), 'r')
    hold on
    plot(t_viz_state, smooth(T_e_viz_state, 10), 'g')
    plot(t_viz_state, smooth(T_w_viz_state, 10), 'b')
    plot(t_viz_state, smooth(T_f_viz_state, 10), 'k')
    set(gca, 'Box', 'off')
    legend('T_S', 'T_E', 'T_W', 'T_F')
    ylabel('Joint Torques (Nm)')
    xlabel('Time (s)')
    ylim([-15 15])
    legend boxoff
    hold off

    % Plot Cost
    figure
    plot(t_viz_state, smooth(Jeffort), 'm')
    set(gca, 'Box', 'off')
    ylabel('Cost (Nm)^2')
    xlabel('Time (s)')
    ylim([0 300])
    hold off

    % Visualize Stroke Trajectories
    lim_y = [-.6 0.1]; lim_x = [-.2 .8]; frames = [130];
    VisualizeDrumStroke
end