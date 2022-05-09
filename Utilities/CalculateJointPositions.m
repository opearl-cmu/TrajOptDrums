% Calculate Joint Positions
XE = L_U.*cos(QS);
YE = L_U.*sin(QS);
XW = L_U.*cos(QS) + L_L.*cos(QE + QS);
YW = L_U.*sin(QS) + L_L.*sin(QE + QS);
XF = L_U.*cos(QS) + L_L.*cos(QE + QS) + L_H.*cos(QE + QS + QW);
YF = L_U.*sin(QS) + L_L.*sin(QE + QS) + L_H.*sin(QE + QS + QW);