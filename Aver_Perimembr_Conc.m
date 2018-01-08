function [Na, Ca] = Aver_Perimembr_Conc(Potentials)
%%
Na_out = 115; % mol/m3 concentration in bulk saline
Na_in = 12 ; % mol/m3 concentration in bulk saline

K_out = 4; % concentration in bulk saline
K_in = 139; % mol/m3 concentration in bulk saline

global Ca_outSR;
Ca_out = Ca_outSR; %1.8  % mol/m3 concentration in bulk saline
Ca_in = 2e-6; % mol/m3 concentration in bulk saline

Cl_out = 120.8; % mol/m3 concentration in bulk saline
Cl_in = 16; % mol/m3 concentration in bulk saline

X_out  = 0; % bid macromolecule with net negative charge
X_in  = 107; % bid macromolecule with net negative charge

Z_Na = 1;  %Valence
Z_K = 1;   %Valence
Z_Ca = 2;  %Valence
Z_Cl = -1; %Valence
Z_X = -1; %Valence
%%
global SR_Temp;
%SR_Temp = 22; %C
T_ = 273.15 + SR_Temp; %K

e_free = 8.85419e-12; %C2 N-1 m-1 or  F?m?1 
R = 8.3144598; % J?mol?1?K?1 Gas constant
Fd = 96485.33289; % C mol?1 Faraday constant
%%
Em_rest = MembrPotentNernst(Na_out, Na_in, K_out, K_in, Ca_out, Ca_in, Cl_out, Cl_in); % trans-membrane potential
V_zero_out = Potentials(1) - Potentials(3); % Outside Zetta potential
V_zero_in = Potentials(2) - Potentials(4) - Em_rest; % Inside Zetta potential

M_out_ = Na_out + K_out + Ca_out + X_out; %approximate out Normality (consoidering alcali metals as interactive species)
e_sol_bulk_T_ = SalPermit_T(T_, M_out_);
Isol_out = 0.5*(Na_out*Z_Na^2 + K_out*Z_K^2 + Ca_out*Z_Ca^2 + Cl_out*Z_Cl^2 + X_out*Z_X^2); % Ionic strength of the Out solution
Isol_in = 0.5*(Na_in*Z_Na^2 + K_in*Z_K^2 + Ca_in*Z_Ca^2 + Cl_in*Z_Cl^2 + X_in*Z_X^2); % Ionic strength of the In solution

%potential distribution%
x = [0.1:0.1:20].*1e-9; % distance from the surface

lamda_out = sqrt(e_free*e_sol_bulk_T_*R*T_/((Fd^2)*2*Isol_out));
V_zero_out_x = exp(-x./lamda_out).*V_zero_out;
Cont_x_out_Na = exp(-V_zero_out_x.*Fd.*Z_Na/(R*T_)).*Na_out;
Cont_x_out_Ca = exp(-V_zero_out_x.*Fd.*Z_Ca/(R*T_)).*Ca_out;

Na = max(Cont_x_out_Na);
Ca = max(Cont_x_out_Ca);

%Na = (max(Cont_x_out_Na) + min(Cont_x_out_Na))/2;
%Ca = (max(Cont_x_out_Ca) + min(Cont_x_out_Ca))/2;

% layerNa = Cont_x_out_Na(find(x<=lamda_out));
% layerCa = Cont_x_out_Ca(find(x<=lamda_out));
% Na = mean(layerNa);
% Ca = mean(layerCa);