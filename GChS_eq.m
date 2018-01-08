function Fout = GChS_eq(Root_in)

F = (1:4); %1:6
VFIout = Root_in(1);
VFIin = Root_in(2);
dVFIout = Root_in(3);
dVFIin = Root_in(4);
%delta_S_out = Root_in(5);
%delta_S_in = Root_in(6);

%--------------------------------------------------------------------------
% temperature 
global SR_Temp;
Tcelsius = SR_Temp; %C
T = 273.15 + Tcelsius; %K
global sigmaPercentKoef

%PARAMETERS
SIGout = -0.016*sigmaPercentKoef;%-0.016; % C*m-1 Intrinsic charge density of the outer bilayer 
SIGin = -0.016/2.7*(1-(sigmaPercentKoef-1)); % C*m-1 Intrinsic charge density of the inner bilayer

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

Vm = MembrPotentNernst(Na_out, Na_in, K_out, K_in, Ca_out, Ca_in, Cl_out, Cl_in);
%Vm = -50*1e-3; %V membrane potential (mV)

%radiiI = [R_Na, R_K, R_Ca,R_Cl];
Radii = RadiiCalcT(T);
R_Na = Radii(1); %0.243e-9; %0.161e-9; % Hydrated ionic radius 
R_K = Radii(2); %0.281e-9; %0.201e-9; % Hydrated ionic radius 
R_Ca = Radii(3); %0.319e-9; %0.168e-9; % Hydrated ionic radius 
R_Cl = Radii(4); %0.246e-9; %0.189e-9; % Hydrated ionic radius 
R_X = 3e-9; % Hydrated ionic radius  % 

e_free = 8.85419e-12; %C2 N-1 m-1 or  F?m?1 
R = 8.3144598; % J?mol?1?K?1 Gas constant
Fd = 96485.33289; % C mol?1 Faraday constant

e_b = 3*e_free; % Permittivity of lipid bilayer (Nanoscale Measurement of the Dielectric Constant of Supported Lipid Bilayers in Aqueous Solutions with Electrostatic Force Microscopy)
delta_b = 3*1e-9; % m  Width of lipid bilayer

M_out = (Na_out + K_out + Ca_out + X_out)/1000; %approximate out Normality (consoidering alcali metals as interactive species)
M_in =  (Na_in + K_in + Ca_in + X_in)/1000; %approximate in Normality (consoidering alcali metals as interactive species)
e_sol_bulk_T = SalPermit_T(T, M_out);
e_sol_bulk = e_sol_bulk_T*e_free; % Permittivity of bulk aqueous solution at temperature T (use the temperature dependence "Dielectric Constant of Saline Water" )

[M_Na, M_Ca] = Aver_Perimembr_Conc(Root_in);
alfa_normal = 1 - 0.255*M_out + 5.151e-2*M_out^2 -6.889e-3*M_out^3;
M_out_Alt = (M_Na + M_Ca + K_out + X_out)/1000; % devision on 1000 done to convert mM into Normality
alfa_Alt = 1 - 0.255*M_out_Alt + 5.151e-2*M_out_Alt^2 -6.889e-3*M_out_Alt^3;
e_S_out = e_sol_bulk*alfa_Alt/alfa_normal; %/10; %  Dielectric constant of outer stern layer
e_S_in = e_sol_bulk*alfa_Alt/alfa_normal; %/10; %  Dielectric constant of inner stern layer
% e_S_out = e_sol_bulk*0.75; %/10; %  Dielectric constant of outer stern layer
% e_S_in = e_sol_bulk*0.75; %/10; %  Dielectric constant of inner stern layer

delta_lipid_out = 0.35e-9; % Hydrated size of the outer polar lipid head groups
delta_lipid_in = 0.35e-9; % Hydrated size of the inner polar lipid head groups
%--------------------------------------------------------------------------

F(1) = (SIGout + (VFIin-VFIout)*e_b/delta_b)^2 - 2*e_sol_bulk*R*T*...
    (Na_out*(exp((-Z_Na*Fd/(R*T))*(VFIout-dVFIout))-1) + ...
    K_out*(exp((-Z_K*Fd/(R*T))*(VFIout-dVFIout))-1) + ...
    Ca_out*(exp((-Z_Ca*Fd/(R*T))*(VFIout-dVFIout))-1) + ...
    Cl_out*(exp((-Z_Cl*Fd/(R*T))*(VFIout-dVFIout))-1) + ...
    X_out*(exp((-Z_X*Fd/(R*T))*(VFIout-dVFIout))-1));

F(2) = (SIGin - (VFIin-VFIout)*e_b/delta_b)^2 - 2*e_sol_bulk*R*T*...
    (Na_in*(exp((-Z_Na*Fd/(R*T))*((VFIin-Vm) - dVFIin))-1) + ...
    K_in*(exp((-Z_K*Fd/(R*T))*((VFIin-Vm) - dVFIin))-1) + ...
    Ca_in*(exp((-Z_Ca*Fd/(R*T))*((VFIin-Vm) - dVFIin))-1) + ...
    Cl_in*(exp((-Z_Cl*Fd/(R*T))*((VFIin-Vm) - dVFIin))-1) + ...
    X_in*(exp((-Z_X*Fd/(R*T))*((VFIin-Vm) - dVFIin))-1));

Aout = (R_Na*Na_out*(exp((-Z_Na*Fd/(R*T))*(VFIout-dVFIout))) + ...
    R_K*K_out*(exp((-Z_K*Fd/(R*T))*(VFIout-dVFIout))) + ...
    R_Ca*Ca_out*(exp((-Z_Ca*Fd/(R*T))*(VFIout-dVFIout))) + ...
    R_Cl*Cl_out*(exp((-Z_Cl*Fd/(R*T))*(VFIout-dVFIout))) + ... 
    R_X*X_out*(exp((-Z_X*Fd/(R*T))*(VFIout-dVFIout))));

Bout = (Na_out*(exp((-Z_Na*Fd/(R*T))*(VFIout-dVFIout))) + ...
    K_out*(exp((-Z_K*Fd/(R*T))*(VFIout-dVFIout))) + ...
    Ca_out*(exp((-Z_Ca*Fd/(R*T))*(VFIout-dVFIout))) + ...
    Cl_out*(exp((-Z_Cl*Fd/(R*T))*(VFIout-dVFIout))) + ... 
    X_out*(exp((-Z_X*Fd/(R*T))*(VFIout-dVFIout))));
%F(5) = delta_S_out - (delta_lipid_out + (Aout)/(Bout));
delta_S_out = (delta_lipid_out + (Aout)/(Bout))

Ain = (R_Na*Na_in*(exp((-Z_Na*Fd/(R*T))*((VFIin-Vm) - dVFIin))) + ...
    R_K*K_in*(exp((-Z_K*Fd/(R*T))*((VFIin-Vm) - dVFIin))) + ...
    R_Ca*Ca_in*(exp((-Z_Ca*Fd/(R*T))*((VFIin-Vm) - dVFIin))) + ...
    R_Cl*Cl_in*(exp((-Z_Cl*Fd/(R*T))*((VFIin-Vm) - dVFIin))) + ...
    R_X*X_in*(exp((-Z_X*Fd/(R*T))*((VFIin-Vm) - dVFIin))));

Bin = (Na_in*(exp((-Z_Na*Fd/(R*T))*((VFIin-Vm) - dVFIin))) + ...
    K_in*(exp((-Z_K*Fd/(R*T))*((VFIin-Vm) - dVFIin))) + ...
    Ca_in*(exp((-Z_Ca*Fd/(R*T))*((VFIin-Vm) - dVFIin))) + ...
    Cl_in*(exp((-Z_Cl*Fd/(R*T))*((VFIin-Vm) - dVFIin))) + ...
    X_in*(exp((-Z_X*Fd/(R*T))*((VFIin-Vm) - dVFIin))));
%F(6) = delta_S_in - (delta_lipid_in + (Ain)/(Bin));
delta_S_in = (delta_lipid_in + (Ain)/(Bin))


F(3) = dVFIout - (delta_S_out/e_S_out)*(SIGout + (VFIin-VFIout)*e_b/delta_b);

F(4) = dVFIin - (delta_S_in/e_S_in)*(SIGin - (VFIin-VFIout)*e_b/delta_b);

Fout = F;

