function EpsBulcSol = SalPermit_T(TTT, M)

TT = TTT - 273.15;
MT=M/1000;

% the temperature dependence "Dielectric Constant of Saline Water"
% Eps_zero = 87.74 - 4.008e-1*TT + 9.398e-4*TT^2 + 1.41e-6*TT^3;
% alfa = 1 - 0.255*MT + 5.151e-2*MT^2 -6.889e-3*MT^3;

% Shuting measurements
Eps_zero = 86.96 - 0.344*TT;
alfa = 1;

EpsBulcSol = Eps_zero*alfa;
