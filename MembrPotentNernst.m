function Em_rest = MembrPotentNernst(Na_out, Na_in, K_out, K_in, Ca_out, Ca_in, Cl_out, Cl_in)

% temperature 
global SR_Temp;
Tcelsius = SR_Temp; %C
T = 273.15 + Tcelsius; %K

Em_rest =  -50*1e-3*(T/298.15); %V membrane potential (mV)