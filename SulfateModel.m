%This model calculates the amount of sodium sulfate required to remove
%barium and strontium in water. Strontium removal is through
%substitution with barium.

tic

%Loading the input data files as mol/L

load('Ba_input')
load('Sr_input')
load('Ca_input')
load('Mg_input')
load('Alk_input')
load('Mg_ends')
load('Ca_ends')
load('Sr_ends')
load('Ba_ends')

%%Defining treatment goals for contaminants (end), solubility constant
%(Ksp for Ba sulfate used), and the Ba:Sr co-precipitation ratio (barium removed
%per strontium

Ca_end = Ca_ends;
Mg_end = Mg_ends;
Ba_end = Ba_ends;
Sr_end = Sr_ends;
s = 5.67;
Ksp = 10^(-9.98);

%Defining the LCA data values
sulfate_AP = 0.00666;
soda_AP = 0.00503;
lime_AP = 0.00153; %Acidification potential (kg SO2 eq per kg soda/lime)
sulfate_EP = 0.00336;
soda_EP = 0.00639;
lime_EP = 0.000597; %Eutrophication potential (kg N eq per kg soda/lime)
sulfate_GWP = 0.648;
soda_GWP = 1.03;
lime_GWP = 0.972; %Global warming potential (kg CO2 eq per kg soda/lime)
sulfate_ODP = 8.83*10^(-8);
soda_ODP = 1.26*10^(-7);
lime_ODP = 7.99*10^(-8); %Ozone depletion potential (kg CFC-11 eq per kg soda/lime)
sulfate_POCP = 0.0552;
soda_POCP = 0.0708;
lime_POCP = 0.0231; %Photochemical ozone creation potential (kg O3 eq per kg soda/lime)
sulfate_PEU = 0.955;
soda_PEU = 1.67;
lime_PEU = 0.717; %Primary energy use (MJ surplus per soda/lime)

%Setting the number of data points in the files (n) and preallocating
%vectors

n = 10000;
Sulfate = zeros(n,1);
ba_add = zeros(n,1);
S_AP = zeros(n,1);
S_EP = zeros(n,1);
S_GWP = zeros(n,1);
S_ODP = zeros(n,1);
S_POCP = zeros(n,1);
S_PEU = zeros(n,1);
Ba_add = zeros(n,1);
Lime = zeros(n,1);
Soda = zeros(n,1);

%Generating a uniform distribution for the "inefficiency factor" found by
%comparing experimental data to model results

sulfate_inefficiency = uni_dist(n, 1.80, 5.77);
lime_inefficiency = uni_dist(n, 1.46, 1.88);
soda_inefficiency = uni_dist(n, 0.77, 1.80);

%% Main loop for Ba and Sr removal

for i = 1:n

%Calculates change in pollutant concentrations

    Ba = Ba_input(i,1) - Ba_end(i,1);    
    if Ba < 0
        Ba = 0;
    end
    Sr = Sr_input(i,1) - Sr_end(i,1);    
    if Sr < 0
        Sr = 0;
    end

%Determines if Ba will limit Sr removal

    if Ba < s * Sr
    
        Ba_add(i,1) = (s * Sr) - Ba;
    
        Ba = Ba + Ba_add(i,1);
      
    else
        
        Ba_add(i,1) = 0;
    
    end    

%Calculates the sulfate required to remove the Ba and Sr

    Sulfate(i,1) = ((Ksp / Ba_end(i,1)) + Sr + Ba) .* sulfate_inefficiency(1,i);
    
%Converts lime and soda ash to kg/m^3

    Sulfate(i,1) = Sulfate(i,1) * 142.04; %Molar mass as sodium sulfate
    
%Calculates emissions impacts per m^3 of water treated
    S_AP(i,1) = (Sulfate(i,1)*sulfate_AP); %Acidification potential
    S_EP(i,1) = (Sulfate(i,1)*sulfate_EP); %Eutrophication potential
    S_GWP(i,1) = (Sulfate(i,1)*sulfate_GWP); %Global warming potential
    S_ODP(i,1) = (Sulfate(i,1)*sulfate_ODP); %Ozone depletion potential
    S_POCP(i,1) = (Sulfate(i,1)*sulfate_POCP); %Photochemical ozone creation potential
    S_PEU(i,1) = (Sulfate(i,1)*sulfate_PEU); %Primary energy use
    
end

%% Main loop for Ca and Mg removal

%Setting solubility product constant for lime soda ash modeling

Ksp = 10^(-8.48);

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

%Calculates carbonate hardness (CCa and CMg) and non-carbonate hardness
%(NCCa and NCMg).

    if 2 * (Ca + Mg) <= Alk_input(i,1)
    
        CCa = Ca;
        NCCa = 0;
        CMg = Mg;
        NCMg = 0;

    elseif 2 * (Ca + Mg) > Alk_input(i,1) && 2 * Ca <= Alk_input(i,1)
    
        CCa = Ca;
        NCCa = 0;
        CMg = (Alk_input(i,1) - 2 * Ca) / 2;
        NCMg = Mg - CMg;
    
    elseif 2 * (Ca + Mg) > Alk_input(i,1) && 2 * Ca > Alk_input(i,1)
    
        CCa = Alk_input(i,1) / 2;
        NCCa = Ca - (0.5 * Alk_input(i,1));
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
    S_AP(i,1) = (Lime(i,1)*lime_AP) + (Soda(i,1)*soda_AP) + S_AP(i,1); %Acidification potential
    S_EP(i,1) = (Lime(i,1)*lime_EP) + (Soda(i,1)*soda_EP) + S_EP(i,1); %Eutrophication potential
    S_GWP(i,1) = (Lime(i,1)*lime_GWP) + (Soda(i,1)*soda_GWP) + S_GWP(i,1); %Global warming potential
    S_ODP(i,1) = (Lime(i,1)*lime_ODP) + (Soda(i,1)*soda_ODP) + S_ODP(i,1); %Ozone depletion potential
    S_POCP(i,1) = (Lime(i,1)*lime_POCP) + (Soda(i,1)*soda_POCP) + S_POCP(i,1); %Photochemical ozone creation potential
    S_PEU(i,1) = (Lime(i,1)*lime_PEU) + (Soda(i,1)*soda_PEU) + S_PEU(i,1); %Primary energy use
    
end

%Saves the sulfate, lime, and soda ash required to a file
csvwrite('Sulfate_Sulfate.csv', Sulfate);
csvwrite('Sulfate_Lime.csv', Lime);
csvwrite('Sulfate_Soda.csv', Soda);
csvwrite('Sulfate_AP.csv', S_AP);
csvwrite('Sulfate_EP.csv', S_EP);
csvwrite('Sulfate_GWP.csv', S_GWP);
csvwrite('Sulfate_ODP.csv', S_ODP);
csvwrite('Sulfate_POCP.csv', S_POCP);
csvwrite('Sulfate_PEU.csv', S_PEU);

%Writes average (median) environmental impact factors to screen
disp(median(S_AP))
disp(median(S_EP))
disp(median(S_GWP))
disp(median(S_ODP))
disp(median(S_POCP))
disp(median(S_PEU))

toc