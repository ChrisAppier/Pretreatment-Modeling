%This model calculates the amount of lime and soda ash required to remove
%calcium, magnesium, barium, and strontium hardness in water. Strontium
%removal is through substitution with calcium.

%Loading the input data files as mol/L

load('Ca_input')
load('Mg_input')
load('Ba_input')
load('Sr_input')
load('Alk_input')
load('Mg_ends')
load('Ca_ends')
load('Sr_ends')
load('Ba_ends')

%Defining treatment goals for contaminants (end), solubility constant
%(Ksp for Ca carbonate used)

Ca_end = Ca_ends;
Mg_end = Mg_ends;
Ba_end = Ba_ends;
Sr_end = Sr_ends;
Ksp = 10^(-8.48);

%Defining the LCA data values
soda_AP = 0.00503;
lime_AP = 0.00153; %Acidification potential (kg SO2 eq per kg soda/lime)
soda_EP = 0.00639;
lime_EP = 0.000597; %Eutrophication potential (kg N eq per kg soda/lime)
soda_GWP = 1.03;
lime_GWP = 0.972; %Global warming potential (kg CO2 eq per kg soda/lime)
soda_ODP = 1.26*10^(-7);
lime_ODP = 7.99*10^(-8); %Ozone depletion potential (kg CFC-11 eq per kg soda/lime)
soda_POCP = 0.0708;
lime_POCP = 0.0231; %Photochemical ozone creation potential (kg O3 eq per kg soda/lime)
soda_PEU = 1.67;
lime_PEU = 0.717; %Primary energy use (MJ surplus per soda/lime)
soda_CAR = 6.43*10^(-08);
lime_CAR = 6.68*10^(-9); %Carcinogenics (CTUh per kg soda/lime)
soda_NCAR = 4.36*10^(-7);
lime_NCAR = 4.30*10^(-8); %Non-carcinogenics (CTUh per kg soda/lime)
soda_RES = 0.001;
lime_RES = 0.000221; %Respiratory effects (kg PM2.5 eq per kg soda/lime)
soda_ETX = 10.9;
lime_ETX = 1.02; %Ecotoxicity (CTUe per kg soda/lime)

%Setting the number of data points in the files (n) and preallocating
%vectors

n = 10000;
Lime = zeros(n:1);
Soda = zeros(n:1);
Sr_end = zeros(n:1);
LSA_AP = zeros(n,1);
LSA_EP = zeros(n,1);
LSA_GWP = zeros(n,1);
LSA_ODP = zeros(n,1);
LSA_POCP = zeros(n,1);
LSA_PEU = zeros(n,1);
LSA_CAR = zeros(n,1);
LSA_NCAR = zeros(n,1);
LSA_RES = zeros(n,1);
LSA_ETX = zeros(n,1);

%Generating a uniform distribution for the "inefficiency factor" found by
%comparing experimental data to model results

lime_inefficiency = uni_dist(n, 1.46, 1.88);
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
    Sr = Sr_input(i,1) * (Ca / Ca_input(i,1)); %Strontium removal is on a 1:1 ratio of percent removal of Ca


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
    
    elseif 2 * (Ca + Ba + Sr + Mg) > Alk_input(i,1) && 2 * (Ca + Ba + Sr) > Alk_input(i,1)
    
        CCa = Alk_input(i,1) / 2;
        NCCa = Ca + Ba + Sr - (0.5 * Alk_input(i,1));
        CMg = 0;
        NCMg = Mg;
    
    end

%Calculates the amount of lime and soda ash required

    Lime(i,1) = CCa + 2*CMg + NCMg + (Ksp / Ca_end(i,1)) .* lime_inefficiency(1,i);

    Soda(i,1) = NCCa + NCMg .* soda_inefficiency(1,i);
    
%Converts lime and soda ash to kg/m^3

    Lime(i,1) = Lime(i,1) * 100.0869;
    
    Soda(i,1) = Soda(i,1) * 105.9888; %molar mass is for anhydrous
    
%Calculates emissions impacts per m^3 of water treated
    LSA_AP(i,1) = (Lime(i,1)*lime_AP) + (Soda(i,1)*soda_AP); %Acidification potential
    LSA_EP(i,1) = (Lime(i,1)*lime_EP) + (Soda(i,1)*soda_EP); %Eutrophication potential
    LSA_GWP(i,1) = (Lime(i,1)*lime_GWP) + (Soda(i,1)*soda_GWP); %Global warming potential
    LSA_ODP(i,1) = (Lime(i,1)*lime_ODP) + (Soda(i,1)*soda_ODP); %Ozone depletion potential
    LSA_POCP(i,1) = (Lime(i,1)*lime_POCP) + (Soda(i,1)*soda_POCP); %Photochemical ozone creation potential
    LSA_PEU(i,1) = (Lime(i,1)*lime_PEU) + (Soda(i,1)*soda_PEU); %Primary energy use
    LSA_CAR(i,1) = (Lime(i,1)*lime_CAR) + (Soda(i,1)*soda_CAR); %Carcinogenics
    LSA_NCAR(i,1) = (Lime(i,1)*lime_NCAR) + (Soda(i,1)*soda_NCAR); %Non-carcinogenics
    LSA_RES(i,1) = (Lime(i,1)*lime_RES) + (Soda(i,1)*soda_RES); %Respiratory effects
    LSA_ETX(i,1) = (Lime(i,1)*lime_ETX) + (Soda(i,1)*soda_ETX); %Ecotoxicity
    
%Calculates the amount of Sr left in water
    
    Sr_end(i,1) = Sr_input(i,1) - Sr;
    if Sr_end(i,1) < 0
        Sr_end(i,1) = 0;
    end
    
end

%% Output

%Saves the lime and soda ash required to a file
csvwrite('LimeSodaAsh_Lime.csv', Lime);
csvwrite('LimeSodaAsh_Soda.csv', Soda);
csvwrite('LimeSodaAsh_AP.csv', LSA_AP);
csvwrite('LimeSodaAsh_EP.csv', LSA_EP);
csvwrite('LimeSodaAsh_GWP.csv', LSA_GWP);
csvwrite('LimeSodaAsh_ODP.csv', LSA_ODP);
csvwrite('LimeSodaAsh_POCP.csv', LSA_POCP);
csvwrite('LimeSodaAsh_PEU.csv', LSA_PEU);
csvwrite('LimeSodaAsh_CAR.csv', LSA_CAR);
csvwrite('LimeSodaAsh_NCAR.csv', LSA_NCAR);
csvwrite('LimeSodaAsh_RES.csv', LSA_RES);
csvwrite('LimeSodaAsh_ETX.csv', LSA_ETX);

%Writes average (median) environmental impact factors to screen
disp(median(LSA_AP))
disp(median(LSA_EP))
disp(median(LSA_GWP))
disp(median(LSA_ODP))
disp(median(LSA_POCP))
disp(median(LSA_PEU))
disp(median(LSA_CAR))
disp(median(LSA_NCAR))
disp(median(LSA_RES))
disp(median(LSA_ETX))