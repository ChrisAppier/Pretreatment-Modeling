%% Limit Calculator

load('SO4_input')

%Creating HCO3 points
HCO3_sample = xlsread('waterdata.xlsx','CG2:CG274');
HCO3_input = lhs_empir(HCO3_sample, n);

%SO4 Based Limits

%Converting SO4_input to mg/l
SO4_input = SO4_input * 96.056 * 1000;

%Constants
mass = SO4_input; %Mass
vol = 1; %Volume
sel = 0.95; %SO4 selectivity
rec = 0.5; %Recovery
Ksp_SO4_Ca = 6.91831*10^(-5); %Ca selectivity
Ksp_SO4_Ba = 1.0842*10^(-10); %Ba selectivity
Ksp_SO4_Sr = 3.44*10^(-7); %Sr selectivity 
Ca_ends = zeros(n,1);
Ba_ends = zeros(n,1);
Sr_ends = zeros(n,1);

%Main loop
for i=1:n
    
    perm = (1-sel)*mass(i);

    rej = (mass(i)-perm)/((1-rec)*vol); %mg/l

    rej = rej / (96*1000); %mol/l

    Ca_sat = Ksp_SO4_Ca/rej; %mol/l
    Ca_ends(i) = Ca_sat*(1-rec)/sel; %mol/l
    
    Ba_sat = Ksp_SO4_Ba/rej; %mol/l
    Ba_ends(i) = Ba_sat*(1-rec)/sel; %mol/l
        
    Sr_sat = Ksp_SO4_Sr/rej; %mol/l
    Sr_ends(i) = Sr_sat*(1-rec)/sel; %mol/l

end

%HCO3 Based Limits, only used if lower than SO4 limits

%Converting HCO3_input to mg/l
HCO3_input = HCO3_input * 1000 * 61; 

%Constants
mass = HCO3_input; %Mass
vol = 1; %Volume
sel = 0.95; %SO4 selectivity
rec = 0.5; %Recovery
Ksp_HCO3_Ca = 3.3*10^(-9); %Ca selectivity
Ksp_HCO3_Ba = 2.58*10^(-9); %Ba selectivity
Ksp_HCO3_Sr = 5.6*10^(-10); %Sr selectivity 

%Main loop
for i=1:n
    
    perm = (1-sel)*mass(i);

    rej = (mass(i)-perm)/((1-rec)*vol); %mg/l

    rej = rej / (61*1000); %mol/l

    Ca_sat = Ksp_HCO3_Ca/rej; %mol/l
    
    if Ca_sat*(1-rec)/sel < Ca_ends(i)
    Ca_ends(i) = Ca_sat*(1-rec)/sel; %mol/l
    end
    
    Ba_sat = Ksp_HCO3_Ba/rej; %mol/l
    
    if Ba_sat*(1-rec)/sel < Ba_ends(i)
    Ba_ends(i) = Ba_sat*(1-rec)/sel; %mol/l
    end
        
    Sr_sat = Ksp_HCO3_Sr/rej; %mol/l
    
    if Sr_sat*(1-rec)/sel < Sr_ends(i)
    Sr_ends(i) = Sr_sat*(1-rec)/sel; %mol/l
    end

end

csvwrite('Ba_ends_test.csv', Ba_ends);
csvwrite('Sr_ends_test.csv', Sr_ends);
csvwrite('Ca_ends_test.csv', Ca_ends);
