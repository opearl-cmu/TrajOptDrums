%-------------------------------+--------------------------+---------------
% Quantity                      | Value                    | Units             
%-------------------------------|--------------------------|---------------
M_T                             =  80;                      % kg               
HEIGHT                          =  1.9;                     % m                
G                               =  9.81;                    % m/s^2            
L_H                             =  0.1668;                  % m                
L_L                             =  0.2924;                  % m                
L_S                             =  0.2080*3/2;              % m                
L_U                             =  0.2554;                  % kg               
M_H                             =  0.0065 * M_T;            % kg               
M_L                             =  0.0187 * M_T;            % kg               
M_S                             =  0.07;                    % kg               
M_U                             =  0.0325 * M_T;            % kg               
I_H                             =  1/12*M_H*L_H^2;          % kg*m^2           
I_L                             =  1/12*M_L*L_L^2;          % kg*m^2           
I_S                             =  1/12*M_S*L_S^2;          % kg*m^2           
I_U                             =  1/12*M_U*L_U^2;          % kg*m^2           
M_ARM                           = M_H + M_L + M_U + M_S;    % k
YD                              = -0.47108;                 % m 