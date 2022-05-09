% Calculate Drum Bead XY-Pose
XB = L_U.*cos(QS) + L_L.*cos(QE+QS) + L_H.*cos(QE+QS+QW) + 0.6666667.*L_S.*cos(QE+QF+QS+QW);
YB = L_U.*sin(QS) + L_L.*sin(QE+QS) + L_H.*sin(QE+QS+QW) + 0.6666667.*L_S.*sin(QE+QF+QS+QW);
YBP = L_U.*cos(QS).*QSP + L_L.*cos(QE+QS).*(QSP+QEP) + L_H.*cos(QE+QS+QW).*(QSP+QEP+QWP) + 0.6666667.*L_S.*cos(QE+QF+QS+QW).*(QSP+QEP+QWP+QFP);

