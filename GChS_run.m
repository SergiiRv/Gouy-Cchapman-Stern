Na_out = 115; % mol/m3 concentration in bulk saline
Na_in = 12 ; % mol/m3 concentration in bulk saline

K_out = 4; % concentration in bulk saline
K_in = 139; % mol/m3 concentration in bulk saline

global Ca_outSR;
Ca_outSR = 1.8;
Ca_out = Ca_outSR; %1.8 % mol/m3 concentration in bulk saline
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

e_free = 8.85419e-12; %C2 N-1 m-1 or  F?m?1 
R = 8.3144598; % J?mol?1?K?1 Gas constant
Fd = 96485.33289; % C mol?1 Faraday constant
%%
VFIout = -25*1e-3; % V (variable)  Potential at the outer bilayer surface
VFIin = -75*1e-3; % V (variable) Potential at the inner bilayer surface
dVFIout = 10*1e-3; % V (variable)  Potential drop at outer bilayer due to stern layer
dVFIin = 10*1e-3; % V (variable)  Potential drop at inner bilayer due to stern layer
delta_S_out = 0.5e-9 ; % m Thickness of stern layer at the outer bilayer
delta_S_in = 0.5e-9; % m  Thickness of stern layer at the inner bilayer

Root_in_ = (1:4); %1:6
Root_in_(1) = VFIout;
Root_in_(2) = VFIin;
Root_in_(3) = dVFIout;
Root_in_(4) = dVFIin;
%Root_in_(5) = delta_S_out;
%Root_in_(6) = delta_S_in;

global SR_Temp;
SR_Temp = 20; %C
T_ = 273.15 + SR_Temp; %K

global sigmaPercentKoef
sigmaPercentKoef = 1%1.25;
%%
disp('the surface potentials: VFIout, VFIin, dVFIout, dVFIin')
fun = @GChS_eq;
options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
Root_solved = fsolve(fun,Root_in_,options)


Em_rest = MembrPotentNernst(Na_out, Na_in, K_out, K_in, Ca_out, Ca_in, Cl_out, Cl_in) % trans-membrane potential
disp('Zetta potential:V_zero_out, V_zero_in')
V_zero_out = Root_solved(1) - Root_solved(3) % Outside Zetta potential
V_zero_in = Root_solved(2) - Root_solved(4) - Em_rest % Inside Zetta potential

M_out_ = Na_out + K_out + Ca_out + X_out; %approximate out Normality (consoidering alcali metals as interactive species)
e_sol_bulk_T_ = SalPermit_T(T_, M_out_);
Isol_out = 0.5*(Na_out*Z_Na^2 + K_out*Z_K^2 + Ca_out*Z_Ca^2 + Cl_out*Z_Cl^2 + X_out*Z_X^2); % Ionic strength of the Out solution
Isol_in = 0.5*(Na_in*Z_Na^2 + K_in*Z_K^2 + Ca_in*Z_Ca^2 + Cl_in*Z_Cl^2 + X_in*Z_X^2); % Ionic strength of the In solution

%potential distribution%
x = [0.05:0.1:15].*1e-9; % distance from the surface

disp('For the temperature'); SR_Temp
lamda_out = sqrt(e_free*e_sol_bulk_T_*R*T_/((Fd^2)*2*Isol_out))
V_zero_out_x = exp(-x./lamda_out).*V_zero_out;
Cont_x_out_Na = exp(-V_zero_out_x.*Fd.*Z_Na/(R*T_)).*Na_out;
Cont_x_out_Ca = exp(-V_zero_out_x.*Fd.*Z_Ca/(R*T_)).*Ca_out;
Max_Na_out = max(Cont_x_out_Na)
Max_Ca_out = max(Cont_x_out_Ca)

% figure('Name','Voltage disribution OUT','NumberTitle','off');
% T_string_short = {strcat('T=', num2str(SR_Temp))};
% plot(x, V_zero_out_x); legend(T_string_short); set(gcf,'color','w'); set(gca,'fontsize',16); 
% xlabel('x (distance), m'); ylabel('Potential V, V');
% figure('Name','Concentration Na disribution OUT','NumberTitle','off');
% plot(x, Cont_x_out_Na);  xlabel('x (distance), m'); ylabel('Concentration, mM'); set(gcf,'color','w'); set(gca,'fontsize',16);
% figure('Name','Concentration Ca disribution OUT','NumberTitle','off');
% plot(x, Cont_x_out_Ca);  xlabel('x (distance), m'); ylabel('Concentration, mM'); set(gcf,'color','w'); set(gca,'fontsize',16);
%%
Temp_begin = 20; % in Celcius
deltaTemp = 1; % in Celcius
Temp_end = 40; % in Celcius
x = [0.05:0.05:20].*1e-9; % distance from the surface
Root_solved_vector = zeros((Temp_end - Temp_begin)/deltaTemp+1,4);
V_zero_vector = zeros((Temp_end - Temp_begin)/deltaTemp+1,2);
Max_Na_Ca_vector = zeros((Temp_end - Temp_begin)/deltaTemp+1,4);
lamda_out_vector = zeros((Temp_end - Temp_begin)/deltaTemp+1,2);
Cont_x_out_Na_vector = zeros((Temp_end - Temp_begin)/deltaTemp+1,length(x));
Cont_x_out_Ca_vector = zeros((Temp_end - Temp_begin)/deltaTemp+1,length(x));
Cont_x_in_Na_vector = zeros((Temp_end - Temp_begin)/deltaTemp+1,length(x));
Cont_x_in_Ca_vector = zeros((Temp_end - Temp_begin)/deltaTemp+1,length(x));
V_zero_out_x_vector = zeros((Temp_end - Temp_begin)/deltaTemp+1,length(x));
V_zero_in_x_vector = zeros((Temp_end - Temp_begin)/deltaTemp+1,length(x));
InMembrane = zeros((Temp_end - Temp_begin)/deltaTemp+1,1); %this is the potential difference between the surfaces of the membrane
InProteine = zeros((Temp_end - Temp_begin)/deltaTemp+1,1); %this is the potential difference between the surfaces of the protein spanning from the surfaceses of the mebrane on 0.5 nm
i = 1;
sigmaPercentKoefVector = [1:0.3/((Temp_end - Temp_begin)/deltaTemp):1.3]; %0.3
for Temp = [Temp_begin:deltaTemp:Temp_end]
    SR_Temp = Temp; %C
    T_ = 273.15 + SR_Temp; %K
    disp('For the temperature'); SR_Temp
    sigmaPercentKoef = sigmaPercentKoefVector(i);
    
    fun = @GChS_eq;
    options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
    Root_solved = fsolve(fun,Root_in_,options);
    Root_solved_vector(i,:) = Root_solved;
    
    Em_rest = MembrPotentNernst(Na_out, Na_in, K_out, K_in, Ca_out, Ca_in, Cl_out, Cl_in); % trans-membrane potential
    disp('Zetta potential:V_zero_out, V_zero_in')
    V_zero_out = Root_solved(1) - Root_solved(3) % Outside Zetta potential
    V_zero_in = Root_solved(2) - Root_solved(4) - Em_rest % Inside Zetta potential
    V_zero_vector(i,:) = [V_zero_out, V_zero_in];
    InMembrane(i,:) = Root_solved(1)-Root_solved(2);

    M_out_ = Na_out + K_out + Ca_out + X_out; %approximate out Normality (consoidering alcali metals as interactive species)
    M_in_ = Na_in + K_in + Ca_in + X_in;
    
    Isol_out = 0.5*(Na_out*Z_Na^2 + K_out*Z_K^2 + Ca_out*Z_Ca^2 + Cl_out*Z_Cl^2 + X_out*Z_X^2); % Ionic strength of the Out solution
    Isol_in = 0.5*(Na_in*Z_Na^2 + K_in*Z_K^2 + Ca_in*Z_Ca^2 + Cl_in*Z_Cl^2 + X_in*Z_X^2); % Ionic strength of the In solution

    %potential distribution%
    
    e_sol_bulk_T_ = SalPermit_T(T_, M_out_);
    lamda_out = sqrt(e_free*e_sol_bulk_T_*R*T_/((Fd^2)*2*Isol_out));
    e_sol_bulk_T_ = SalPermit_T(T_, M_in_);
    lamda_in = sqrt(e_free*e_sol_bulk_T_*R*T_/((Fd^2)*2*Isol_in));
    V_zero_out_x = exp(-x./lamda_out).*V_zero_out;
    V_zero_in_x = exp(-x./lamda_in).*V_zero_in;
    InProteine(i,:) = V_zero_out_x(x==0.1*1.0e-9) - (V_zero_in_x(x==0.1*1.0e-9)+Em_rest);
    Cont_x_out_Na = exp(-V_zero_out_x.*Fd.*Z_Na/(R*T_)).*Na_out;
    Cont_x_out_Ca = exp(-V_zero_out_x.*Fd.*Z_Ca/(R*T_)).*Ca_out;
    Cont_x_in_Na = exp(-V_zero_in_x.*Fd.*Z_Na/(R*T_)).*Na_in;
    Cont_x_in_Ca = exp(-V_zero_in_x.*Fd.*Z_Ca/(R*T_)).*Ca_in;
    Max_Na_out = max(Cont_x_out_Na);
    Max_Ca_out = max(Cont_x_out_Ca);
    Max_Na_in = max(Cont_x_in_Na);
    Max_Ca_in = max(Cont_x_in_Ca);
    
    lamda_out_vector(i,:) = [lamda_out, lamda_in];
    Max_Na_Ca_vector(i,:) = [Max_Na_out, Max_Ca_out, Max_Na_in, Max_Ca_in];
    Cont_x_out_Na_vector(i,:) = Cont_x_out_Na;
    Cont_x_out_Ca_vector(i,:) = Cont_x_out_Ca;
    Cont_x_in_Na_vector(i,:) = Cont_x_in_Na;
    Cont_x_in_Ca_vector(i,:) = Cont_x_in_Ca;
    V_zero_out_x_vector(i,:) = V_zero_out_x;
    V_zero_in_x_vector(i,:) = V_zero_in_x;
    %----------------------------------------------------------------------
    i = i+1;
end
Temp = [Temp_begin:deltaTemp:Temp_end];
T_string = cell(1, length(Temp));
for ii = 1:length(Temp)
    T_string(ii) = {strcat('T=', num2str(Temp(ii)))};
end


figure('Name','Voltage disribution OUT: for different temp','NumberTitle','off');
plot(x, V_zero_out_x_vector); legend(T_string); xlabel('x (distance), m'); ylabel('Potential V, V'); set(gcf,'color','w'); set(gca,'fontsize',14);
saveas(gcf,'Voltage disribution OUT for different temp') 
    figure('Name','Voltage in Membrane and protein: for different temp','NumberTitle','off');
    plot(Temp, InMembrane, Temp, InProteine); xlabel('temperature, C'); ylabel('DeltaPotential V, V'); set(gcf,'color','w'); set(gca,'fontsize',14); legend({'InMembrane', 'InProteine'});
    saveas(gcf,'Voltage in Membrane and protein for different temp') 
figure('Name','Concentration Na disribution OUT: for different temp','NumberTitle','off');
plot(x, Cont_x_out_Na_vector); legend(T_string); xlabel('x (distance), m'); ylabel('Concentration, mM'); set(gcf,'color','w'); set(gca,'fontsize',14);
saveas(gcf,'Concentration Na disribution OUT for different temp')
figure('Name','Concentration Ca disribution OUT: for different temp','NumberTitle','off');
plot(x, Cont_x_out_Ca_vector); legend(T_string); xlabel('x (distance), m'); ylabel('Concentration, mM'); set(gcf,'color','w'); set(gca,'fontsize',14);
saveas(gcf,'Concentration Ca disribution OUT for different temp')
figure('Name','Lambda for different temp','NumberTitle','off');
plot(Temp, lamda_out_vector); xlabel('temperature, C'); ylabel('Lambda, m'); set(gcf,'color','w'); set(gca,'fontsize',14);
saveas(gcf,'Lambda for different temp')
figure('Name','Max_Na Concentration for different temp','NumberTitle','off');
plot(Temp, Max_Na_Ca_vector(:,1)); xlabel('temperature, C'); ylabel('Concentration, mM'); set(gcf,'color','w'); set(gca,'fontsize',14);
saveas(gcf,'Max_Na Concentration for different temp')
figure('Name','Max_Ca Concentration for different temp','NumberTitle','off');
plot(Temp, Max_Na_Ca_vector(:,2)); xlabel('temperature, C'); ylabel('Concentration, mM'); set(gcf,'color','w'); set(gca,'fontsize',14);
saveas(gcf,'Max_Ca Concentration for different temp')
figure('Name','V_zero_Out for different temp','NumberTitle','off');
plot(Temp, V_zero_vector(:,1)); xlabel('temperature, C'); ylabel('Max Surface Potential (Zetta potential), V'); set(gcf,'color','w'); set(gca,'fontsize',14);
saveas(gcf,'V_zero_Out for different temp')
figure('Name','V_zero_In for different temp','NumberTitle','off');
plot(Temp, V_zero_vector(:,2)); xlabel('temperature, C'); ylabel('Max Surface Potential (Zetta potential), V'); set(gcf,'color','w'); set(gca,'fontsize',14);
saveas(gcf,'V_zero_In for different temp')

filename = 'C:\MathLab_test_folder\Gouy-Chapman-Stern theory\Run for Retzius paper\GChS_result_T20to60_VAR_Sig1_0to1_3_detailed.xlsx';
A = {'Temperature'; 'x  distance'}; sheet = 'V_zero_out_x'; xlRange = 'A1'; xlswrite(filename,A,sheet,xlRange);
A = Temp; sheet = 'V_zero_out_x'; xlRange = 'B1'; xlswrite(filename,A,sheet,xlRange);
A = x'; sheet = 'V_zero_out_x'; xlRange = 'A3'; xlswrite(filename,A,sheet,xlRange);
A = V_zero_out_x_vector'; sheet = 'V_zero_out_x'; xlRange = 'B3'; xlswrite(filename,A,sheet,xlRange);

A = {'Temperature'; 'x  distance'}; sheet = 'V_zero_in_x'; xlRange = 'A1'; xlswrite(filename,A,sheet,xlRange);
A = Temp; sheet = 'V_zero_in_x'; xlRange = 'B1'; xlswrite(filename,A,sheet,xlRange);
A = x'; sheet = 'V_zero_in_x'; xlRange = 'A3'; xlswrite(filename,A,sheet,xlRange);
A = V_zero_in_x_vector'; sheet = 'V_zero_in_x'; xlRange = 'B3'; xlswrite(filename,A,sheet,xlRange);

A = {'Temperature'; 'x  distance'}; sheet = 'Cont_x_out_Na'; xlRange = 'A1'; xlswrite(filename,A,sheet,xlRange);
A = Temp; sheet = 'Cont_x_out_Na'; xlRange = 'B1'; xlswrite(filename,A,sheet,xlRange);
A = x'; sheet = 'Cont_x_out_Na'; xlRange = 'A3'; xlswrite(filename,A,sheet,xlRange);
A = Cont_x_out_Na_vector'; sheet = 'Cont_x_out_Na'; xlRange = 'B3'; xlswrite(filename,A,sheet,xlRange);

A = {'Temperature'; 'x  distance'}; sheet = 'Cont_x_out_Ca'; xlRange = 'A1'; xlswrite(filename,A,sheet,xlRange);
A = Temp; sheet = 'Cont_x_out_Ca'; xlRange = 'B1'; xlswrite(filename,A,sheet,xlRange);
A = x'; sheet = 'Cont_x_out_Ca'; xlRange = 'A3'; xlswrite(filename,A,sheet,xlRange);
A = Cont_x_out_Ca_vector'; sheet = 'Cont_x_out_Ca'; xlRange = 'B3'; xlswrite(filename,A,sheet,xlRange);

A = {'Temperature'}; sheet = 'lamda_out'; xlRange = 'A1'; xlswrite(filename,A,sheet,xlRange);
A = Temp'; sheet = 'lamda_out'; xlRange = 'A2'; xlswrite(filename,A,sheet,xlRange);
A = {'lamda_out', 'lamda_in'}; sheet = 'lamda_out'; xlRange = 'B1'; xlswrite(filename,A,sheet,xlRange);
A = lamda_out_vector; sheet = 'lamda_out'; xlRange = 'B2'; xlswrite(filename,A,sheet,xlRange);

A = {'Temperature'}; sheet = 'InMembrane-Proteine'; xlRange = 'A1'; xlswrite(filename,A,sheet,xlRange);
A = Temp'; sheet = 'InMembrane-Proteine'; xlRange = 'A2'; xlswrite(filename,A,sheet,xlRange);
A = {'InMembrane', 'InProteine'}; sheet = 'InMembrane-Proteine'; xlRange = 'B1'; xlswrite(filename,A,sheet,xlRange);
A = [InMembrane, InProteine]; sheet = 'InMembrane-Proteine'; xlRange = 'B2'; xlswrite(filename,A,sheet,xlRange);


A = {'Temperature', 'Max_Na_Out', 'Max_Ca_Out', 'Max_Na_In', 'Max_Ca_In'}; sheet = 'Max_Na_Ca_Out'; xlRange = 'A1'; xlswrite(filename,A,sheet,xlRange);
A = Temp'; sheet = 'Max_Na_Ca_Out'; xlRange = 'A2'; xlswrite(filename,A,sheet,xlRange);
A = Max_Na_Ca_vector; sheet = 'Max_Na_Ca_Out'; xlRange = 'B2'; xlswrite(filename,A,sheet,xlRange);


A = {'Temperature', 'V_zero_out', 'V_zero_in'}; sheet = 'Max_V_zero_Out-In'; xlRange = 'A1'; xlswrite(filename,A,sheet,xlRange);
A = Temp'; sheet = 'Max_V_zero_Out-In'; xlRange = 'A2'; xlswrite(filename,A,sheet,xlRange);
A = V_zero_vector; sheet = 'Max_V_zero_Out-In'; xlRange = 'B2'; xlswrite(filename,A,sheet,xlRange);