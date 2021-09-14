%This model calculates the amount of caustic soda required to remove
%calcium, magnesium, barium, and strontium hardness in water.

%Loading the input data files as mol/L

load('Ca_input')
load('Mg_input')
load('Ba_input')
load('Sr_input')
load('Alk_input')
load('inefficiency')
load('Mg_ends')
load('Ca_ends')
load('Sr_ends')
load('Ba_ends')

%Defining treatment goals for contaminants (mol/L), solubility constant
%(Ksp for Ca carbonate used), and the kWh energy required per mol/L of caustic soda (e_fac).

Ca_end = Ca_ends;
Mg_end = Mg_ends;
Ba_end = Ba_ends;
Sr_end = Sr_ends;
Ksp = 10^(-8.48);
e_fac  = 0.159988;

%Defining the LCA data values
soda_AP = 0.00503;
energy_AP = 0.000564; %Acidification potential (kg SO2 eq per kg caustic soda or kWh energy used)
soda_EP = 0.00639;
energy_EP = 0.00148; %Eutrophication potential (kg N eq per kg caustic soda or kWh energy used)
soda_GWP = 1.03;
energy_GWP = 0.182; %Global warming potential (kg CO2 eq per kg caustic soda or kWh energy used)
soda_ODP = 1.26*10^(-7);
energy_ODP = 1.51*10^(-8); %Ozone depletion potential (kg CFC-11 eq per kg caustic soda or kWh energy used)
soda_POCP = 0.0708;
energy_POCP = 0.00483; %Photochemical ozone creation potential (kg O3 eq per kg caustic soda or kWh energy used)
soda_PEU = 1.67;
energy_PEU = 0.131; %Primary energy use (MJ surplus per kg caustic soda or kWh energy used)
soda_CAR = 6.43*10^(-8);
energy_CAR = 1.55*10^(-8); %Carcinogenics (CTUh per kg caustic soda or kWh energy used)
soda_NCAR = 4.36*10^(-7);
energy_NCAR = 5.22*10^(-8); %Non-carcinogenics (CTUh per kg caustic soda or kWh energy used)
soda_RES = 0.000496;
energy_RES = 0.000556; %Respiratory effects (kg PM2.5 eq per kg caustic soda or kWh energy used)
soda_ETX = 10.9;
energy_ETX = 1.46; %Ecotoxicity (CTUe per kg caustic soda or kWh energy used)
    
%Setting the number of data points in the files (n) and preallocating
%vectors

n = 10000;
Caustic = zeros(n,1);
Soda = zeros(n,1);
Energy = zeros(n,1);
CS_AP = zeros(n,1);
CS_EP = zeros(n,1);
CS_GWP = zeros(n,1);
CS_ODP = zeros(n,1);
CS_POCP = zeros(n,1);
CS_PEU = zeros(n,1);
CS_CAR = zeros(n,1);
CS_NCAR = zeros(n,1);
CS_RES = zeros(n,1);
CS_ETX = zeros(n,1);

%Generating a uniform distribution for the "inefficiency factor" found by
%comparing experimental data to model results

caustic_inefficiency = uni_dist(n, 1.64, 2.13);
soda_inefficiency = uni_dist(n, 0.77, 1.80);

%% Main loop

for i = 1:n

%Calculates change in pollutant concentrations

    Ca = Ca_input(i,1) - Ca_end(i,1);
    if Ca < 0
        Ca = 0;
    end
    Mg = Mg_input(i,1) - Mg_end(i,1);    
    if Mg < 0
        Mg = 0;
    end
    Ba = Ba_input(i,1) - Ba_end(i,1);    
    if Ba < 0
        Ba = 0;
    end
    Sr = Sr_input(i,1) - Sr_end(i,1);    
    if Sr < 0
        Sr = 0;
    end

%Calculates carbonate hardness (CCa and CMg) and non-carbonate hardness
%(NCCa and NCMg). Ba and Sr behave similar to Ca and are all combined under
%the CCa and NCCa variables

    if 2 * (Ca + Ba + Sr + Mg) <= Alk_input(i,1)
    
        CCa = Ca + Ba + Sr;
        NCCa = 0;
        CMg = Mg;
        NCMg = 0;

    elseif 2 * (Ca + Ba + Sr + Mg) > Alk_input(i,1) && 2 * (Ca + Ba + Sr) <= Alk_input(i,1)
    
        CCa = Ca + Ba + Sr;
        NCCa = 0;
        CMg = (Alk_input(i,1) - 2 * (Ca + Ba + Sr)) / 2;
        NCMg = Mg - CMg;
    
    elseif 2* (Ca + Ba + Sr + Mg) > Alk_input(i,1) && 2* (Ca + Ba + Sr) > Alk_input(i,1)
    
        CCa = Alk_input(i,1) / 2;
        NCCa = Ca + Ba + Sr - (0.5 * Alk_input(i,1));
        CMg = 0;
        NCMg = Mg;
    
    end
    
%Calculates the amount of caustic soda, soda ash required, and energy
%required

    Caustic(i,1) = (2*CCa + 4*CMg + 2*NCMg + (Ksp / Ca_end(i))) .* caustic_inefficiency(1,i);

    Soda(i,1) = NCCa .* soda_inefficiency(1,i);
    
    Energy(i,1) = Caustic(i,1) * e_fac;
    
%Converts caustic soda and soda ash to kg/m^3

    Caustic(i,1) = Caustic(i,1) * 39.9971;
    
    Soda(i,1) = Soda(i,1) * 105.9888; %molar mass is for anhydrous
    
%Calculates emissions impacts per m^3 of water treated
    CS_AP(i,1) = (Soda(i,1)*soda_AP) + (Energy(i,1)*energy_AP); %Acidification potential
    CS_EP(i,1) = (Soda(i,1)*soda_EP) + (Energy(i,1)*energy_EP); %Eutrophication potential
    CS_GWP(i,1) = (Soda(i,1)*soda_GWP) + (Energy(i,1)*energy_GWP); %Global warming potential
    CS_ODP(i,1) = (Soda(i,1)*soda_ODP) + (Energy(i,1)*energy_ODP); %Ozone depletion potential
    CS_POCP(i,1) = (Soda(i,1)*soda_POCP) + (Energy(i,1)*energy_POCP); %Photochemical ozone creation potential
    CS_PEU(i,1) = (Soda(i,1)*soda_PEU) + (Energy(i,1)*energy_PEU); %Fossil Fuel Depletion
    CS_CAR(i,1) = (Soda(i,1)*soda_CAR) + (Energy(i,1)*energy_CAR); %Carcinogenics
    CS_NCAR(i,1) = (Soda(i,1)*soda_NCAR) + (Energy(i,1)*energy_NCAR); %Non-carcinogenics
    CS_RES(i,1) = (Soda(i,1)*soda_RES) + (Energy(i,1)*energy_RES); %Respiratory effects
    CS_ETX(i,1) = (Soda(i,1)*soda_ETX) + (Energy(i,1)*energy_ETX); %Ecotoxicity

end

%% Output

%Saves the energy, caustic soda, and soda ash required to a file ***need to
%add environmental impact factors to an excel file***
csvwrite('CausticSoda_Caustic.csv', Caustic);
csvwrite('CausticSoda_Soda.csv', Soda);
csvwrite('CausticSoda_Energy.csv', Energy);
csvwrite('CausticSoda_AP.csv', CS_AP);
csvwrite('CausticSoda_EP.csv', CS_EP);
csvwrite('CausticSoda_GWP.csv', CS_GWP);
csvwrite('CausticSoda_ODP.csv', CS_ODP);
csvwrite('CausticSoda_POCP.csv', CS_POCP);
csvwrite('CausticSoda_PEU.csv', CS_PEU);
csvwrite('CausticSoda_CAR.csv', CS_CAR);
csvwrite('CausticSoda_NCAR.csv', CS_NCAR);
csvwrite('CausticSoda_RES.csv', CS_RES);
csvwrite('CausticSoda_ETX.csv', CS_ETX);
csvwrite('Mg.csv', Mg);

%Writes average (median) environmental impact factors to screen
disp(median(CS_AP))
disp(median(CS_EP))
disp(median(CS_GWP))
disp(median(CS_ODP))
disp(median(CS_POCP))
disp(median(CS_PEU))
disp(median(CS_CAR))
disp(median(CS_NCAR))
disp(median(CS_RES))
disp(median(CS_ETX))