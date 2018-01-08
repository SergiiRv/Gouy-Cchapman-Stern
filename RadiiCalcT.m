function radiiI = RadiiCalcT(TT)

radiiI = [1:4];
%TT is a copy of T but have only local properties 
k = 1.38e-23; %J/K

R_Na_0 = 0.091e-9; %m
R_K_0 = 0.143e-9; %m
R_Cl_0 = 0.195e-9; %m
R_Ca_0 = 0.085e-9; %m

dR_Na = 0.152e-9/2; %m
dR_K = 0.138e-9/1.8; %m
dR_Cl = 0.124e-9/2; %m
dR_Ca = 0.161e-9/2.2; %m

A_Na = 0.0397;
A_K = 0.0114;
A_Cl = 0.0049;
A_Ca = 2.37;

U_Na = -2.23e-20; %J
U_K = -2.53e-20; %J
U_Cl = -2.82e-20; %J
U_Ca = -0.6e-20; %J

R_Na = R_Na_0 + dR_Na*(A_Na*exp(-U_Na/(k*TT)))^(1/3);
R_K = R_K_0 + dR_K*(A_K*exp(-U_K/(k*TT)))^(1/3);
R_Ca = R_Ca_0 + dR_Ca*(A_Ca*exp(-U_Ca/(k*TT)))^(1/3);
R_Cl = R_Cl_0 + dR_Cl*(A_Cl*exp(-U_Cl/(k*TT)))^(1/3);
radiiI = [R_Na, R_K, R_Ca, R_Cl];
