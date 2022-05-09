% CALCULATE COST GRADIENTS
DJDT = zeros(size(t));
DJDQS = zeros(size(t));
DJDQE = zeros(size(t));
DJDQW = zeros(size(t));
DJDQF = zeros(size(t));
DJDQSP = zeros(size(t));
DJDQEP = zeros(size(t));
DJDQWP = zeros(size(t));
DJDQFP = zeros(size(t));
DJDT_S = 2.*T_S;
DJDT_E = 2.*T_E;
DJDT_W = 2.*T_W;
DJDT_F = 2.*T_F;
DJDF_D = zeros(size(t));





