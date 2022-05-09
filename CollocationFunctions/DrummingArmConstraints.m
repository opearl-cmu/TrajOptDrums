function [Cp_in, Cp_eq, Cp_in_Grad, Cp_eq_Grad] = DrummingArmConstraints(t, z, u, trackXB)
% Load Model Parameters
SetParams

% Unpack State
QS = z(1,:);
QE = z(2,:);
QW = z(3,:);
QF = z(4,:);
QSP = z(5,:);
QEP = z(6,:);
QWP = z(7,:);
QFP = z(8,:);

% Unpack Force
F_D = u(5,:);

% Calculate YB from Model State
CalculateDrumBeadPose

% Calculate Joint Centers
CalculateJointPositions

% Prescribe YB, XB, and F_D
XB_track = Prescribe_XB(t); XB_track_err = 0.01; % m
YB_track = Prescribe_YB(t); YB_track_err = 0.01; % m
F_D_track = Prescribe_F_D(t); F_D_track_err = 0.1; % Newtons

% Other misc. constraint settings/bounds
offsetYDlim = -0.003; XB_low = 0.4; XB_upp = 0.65;

% Set Equality Constraints and Gradients (don't use equality constraints in optitraj toolbox because IPOPT takes everything as inequality constraint inputs)
Cp_eq = [];
Cp_eq_Grad = [];

% Set Inequality Constraints & Gradients
if trackXB == true
    Cp_in = [YD+offsetYDlim-YE; YD+offsetYDlim-YW; YD+offsetYDlim-YF; XB-XB_upp; XB_low-XB; YB-(YB_track+YB_track_err); (YB_track-YB_track_err)-YB; F_D-(F_D_track+F_D_track_err); (F_D_track-F_D_track_err)-F_D;   XB-(XB_track+XB_track_err); (XB_track-XB_track_err)-XB];
    CalculateInequalityGradients_TrackXB;
    Cp_in_Grad = permute( cat(3, [DCP_IN_ONEDT; DCP_IN_ONEDQS; DCP_IN_ONEDQE; DCP_IN_ONEDQW; DCP_IN_ONEDQF; DCP_IN_ONEDQSP; DCP_IN_ONEDQEP; DCP_IN_ONEDQWP; DCP_IN_ONEDQFP; DCP_IN_ONEDT_S; DCP_IN_ONEDT_E; DCP_IN_ONEDT_W; DCP_IN_ONEDT_F; DCP_IN_ONEDF_D], [DCP_IN_TWODT; DCP_IN_TWODQS; DCP_IN_TWODQE; DCP_IN_TWODQW; DCP_IN_TWODQF; DCP_IN_TWODQSP; DCP_IN_TWODQEP; DCP_IN_TWODQWP; DCP_IN_TWODQFP; DCP_IN_TWODT_S; DCP_IN_TWODT_E; DCP_IN_TWODT_W; DCP_IN_TWODT_F; DCP_IN_TWODF_D], [DCP_IN_THREEDT; DCP_IN_THREEDQS; DCP_IN_THREEDQE; DCP_IN_THREEDQW; DCP_IN_THREEDQF; DCP_IN_THREEDQSP; DCP_IN_THREEDQEP; DCP_IN_THREEDQWP; DCP_IN_THREEDQFP; DCP_IN_THREEDT_S; DCP_IN_THREEDT_E; DCP_IN_THREEDT_W; DCP_IN_THREEDT_F; DCP_IN_THREEDF_D], [DCP_IN_FOURDT; DCP_IN_FOURDQS; DCP_IN_FOURDQE; DCP_IN_FOURDQW; DCP_IN_FOURDQF; DCP_IN_FOURDQSP; DCP_IN_FOURDQEP; DCP_IN_FOURDQWP; DCP_IN_FOURDQFP; DCP_IN_FOURDT_S; DCP_IN_FOURDT_E; DCP_IN_FOURDT_W; DCP_IN_FOURDT_F; DCP_IN_FOURDF_D], [DCP_IN_FIVEDT; DCP_IN_FIVEDQS; DCP_IN_FIVEDQE; DCP_IN_FIVEDQW; DCP_IN_FIVEDQF; DCP_IN_FIVEDQSP; DCP_IN_FIVEDQEP; DCP_IN_FIVEDQWP; DCP_IN_FIVEDQFP; DCP_IN_FIVEDT_S; DCP_IN_FIVEDT_E; DCP_IN_FIVEDT_W; DCP_IN_FIVEDT_F; DCP_IN_FIVEDF_D], [DCP_IN_SIXDT; DCP_IN_SIXDQS; DCP_IN_SIXDQE; DCP_IN_SIXDQW; DCP_IN_SIXDQF; DCP_IN_SIXDQSP; DCP_IN_SIXDQEP; DCP_IN_SIXDQWP; DCP_IN_SIXDQFP; DCP_IN_SIXDT_S; DCP_IN_SIXDT_E; DCP_IN_SIXDT_W; DCP_IN_SIXDT_F; DCP_IN_SIXDF_D], [DCP_IN_SEVENDT; DCP_IN_SEVENDQS; DCP_IN_SEVENDQE; DCP_IN_SEVENDQW; DCP_IN_SEVENDQF; DCP_IN_SEVENDQSP; DCP_IN_SEVENDQEP; DCP_IN_SEVENDQWP; DCP_IN_SEVENDQFP; DCP_IN_SEVENDT_S; DCP_IN_SEVENDT_E; DCP_IN_SEVENDT_W; DCP_IN_SEVENDT_F; DCP_IN_SEVENDF_D], [DCP_IN_EIGHTDT; DCP_IN_EIGHTDQS; DCP_IN_EIGHTDQE; DCP_IN_EIGHTDQW; DCP_IN_EIGHTDQF; DCP_IN_EIGHTDQSP; DCP_IN_EIGHTDQEP; DCP_IN_EIGHTDQWP; DCP_IN_EIGHTDQFP; DCP_IN_EIGHTDT_S; DCP_IN_EIGHTDT_E; DCP_IN_EIGHTDT_W; DCP_IN_EIGHTDT_F; DCP_IN_EIGHTDF_D], [DCP_IN_NINEDT; DCP_IN_NINEDQS; DCP_IN_NINEDQE; DCP_IN_NINEDQW; DCP_IN_NINEDQF; DCP_IN_NINEDQSP; DCP_IN_NINEDQEP; DCP_IN_NINEDQWP; DCP_IN_NINEDQFP; DCP_IN_NINEDT_S; DCP_IN_NINEDT_E; DCP_IN_NINEDT_W; DCP_IN_NINEDT_F; DCP_IN_NINEDF_D], [DCP_IN_TENDT; DCP_IN_TENDQS; DCP_IN_TENDQE; DCP_IN_TENDQW; DCP_IN_TENDQF; DCP_IN_TENDQSP; DCP_IN_TENDQEP; DCP_IN_TENDQWP; DCP_IN_TENDQFP; DCP_IN_TENDT_S; DCP_IN_TENDT_E; DCP_IN_TENDT_W; DCP_IN_TENDT_F; DCP_IN_TENDF_D], [DCP_IN_ELEVENDT; DCP_IN_ELEVENDQS; DCP_IN_ELEVENDQE; DCP_IN_ELEVENDQW; DCP_IN_ELEVENDQF; DCP_IN_ELEVENDQSP; DCP_IN_ELEVENDQEP; DCP_IN_ELEVENDQWP; DCP_IN_ELEVENDQFP; DCP_IN_ELEVENDT_S; DCP_IN_ELEVENDT_E; DCP_IN_ELEVENDT_W; DCP_IN_ELEVENDT_F; DCP_IN_ELEVENDF_D]), [3, 1, 2]);
else
    Cp_in = [YD+offsetYDlim-YE; YD+offsetYDlim-YW; YD+offsetYDlim-YF; XB-XB_upp; XB_low-XB; YB-(YB_track+YB_track_err); (YB_track-YB_track_err)-YB; F_D-(F_D_track+F_D_track_err); (F_D_track-F_D_track_err)-F_D];
    CalculateInequalityGradients;
    Cp_in_Grad = permute( cat(3, [DCP_IN_ONEDT; DCP_IN_ONEDQS; DCP_IN_ONEDQE; DCP_IN_ONEDQW; DCP_IN_ONEDQF; DCP_IN_ONEDQSP; DCP_IN_ONEDQEP; DCP_IN_ONEDQWP; DCP_IN_ONEDQFP; DCP_IN_ONEDT_S; DCP_IN_ONEDT_E; DCP_IN_ONEDT_W; DCP_IN_ONEDT_F; DCP_IN_ONEDF_D], [DCP_IN_TWODT; DCP_IN_TWODQS; DCP_IN_TWODQE; DCP_IN_TWODQW; DCP_IN_TWODQF; DCP_IN_TWODQSP; DCP_IN_TWODQEP; DCP_IN_TWODQWP; DCP_IN_TWODQFP; DCP_IN_TWODT_S; DCP_IN_TWODT_E; DCP_IN_TWODT_W; DCP_IN_TWODT_F; DCP_IN_TWODF_D], [DCP_IN_THREEDT; DCP_IN_THREEDQS; DCP_IN_THREEDQE; DCP_IN_THREEDQW; DCP_IN_THREEDQF; DCP_IN_THREEDQSP; DCP_IN_THREEDQEP; DCP_IN_THREEDQWP; DCP_IN_THREEDQFP; DCP_IN_THREEDT_S; DCP_IN_THREEDT_E; DCP_IN_THREEDT_W; DCP_IN_THREEDT_F; DCP_IN_THREEDF_D], [DCP_IN_FOURDT; DCP_IN_FOURDQS; DCP_IN_FOURDQE; DCP_IN_FOURDQW; DCP_IN_FOURDQF; DCP_IN_FOURDQSP; DCP_IN_FOURDQEP; DCP_IN_FOURDQWP; DCP_IN_FOURDQFP; DCP_IN_FOURDT_S; DCP_IN_FOURDT_E; DCP_IN_FOURDT_W; DCP_IN_FOURDT_F; DCP_IN_FOURDF_D], [DCP_IN_FIVEDT; DCP_IN_FIVEDQS; DCP_IN_FIVEDQE; DCP_IN_FIVEDQW; DCP_IN_FIVEDQF; DCP_IN_FIVEDQSP; DCP_IN_FIVEDQEP; DCP_IN_FIVEDQWP; DCP_IN_FIVEDQFP; DCP_IN_FIVEDT_S; DCP_IN_FIVEDT_E; DCP_IN_FIVEDT_W; DCP_IN_FIVEDT_F; DCP_IN_FIVEDF_D], [DCP_IN_SIXDT; DCP_IN_SIXDQS; DCP_IN_SIXDQE; DCP_IN_SIXDQW; DCP_IN_SIXDQF; DCP_IN_SIXDQSP; DCP_IN_SIXDQEP; DCP_IN_SIXDQWP; DCP_IN_SIXDQFP; DCP_IN_SIXDT_S; DCP_IN_SIXDT_E; DCP_IN_SIXDT_W; DCP_IN_SIXDT_F; DCP_IN_SIXDF_D], [DCP_IN_SEVENDT; DCP_IN_SEVENDQS; DCP_IN_SEVENDQE; DCP_IN_SEVENDQW; DCP_IN_SEVENDQF; DCP_IN_SEVENDQSP; DCP_IN_SEVENDQEP; DCP_IN_SEVENDQWP; DCP_IN_SEVENDQFP; DCP_IN_SEVENDT_S; DCP_IN_SEVENDT_E; DCP_IN_SEVENDT_W; DCP_IN_SEVENDT_F; DCP_IN_SEVENDF_D], [DCP_IN_EIGHTDT; DCP_IN_EIGHTDQS; DCP_IN_EIGHTDQE; DCP_IN_EIGHTDQW; DCP_IN_EIGHTDQF; DCP_IN_EIGHTDQSP; DCP_IN_EIGHTDQEP; DCP_IN_EIGHTDQWP; DCP_IN_EIGHTDQFP; DCP_IN_EIGHTDT_S; DCP_IN_EIGHTDT_E; DCP_IN_EIGHTDT_W; DCP_IN_EIGHTDT_F; DCP_IN_EIGHTDF_D], [DCP_IN_NINEDT; DCP_IN_NINEDQS; DCP_IN_NINEDQE; DCP_IN_NINEDQW; DCP_IN_NINEDQF; DCP_IN_NINEDQSP; DCP_IN_NINEDQEP; DCP_IN_NINEDQWP; DCP_IN_NINEDQFP; DCP_IN_NINEDT_S; DCP_IN_NINEDT_E; DCP_IN_NINEDT_W; DCP_IN_NINEDT_F; DCP_IN_NINEDF_D]), [3, 1, 2]);
end

end