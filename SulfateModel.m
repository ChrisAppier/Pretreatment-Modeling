%This model calculates the amount of sodium sulfate required to remove
%barium and strontium in water. Strontium removal is through
%substitution with barium.

clc 
clear
tic

%Loading the input data files as mol/L

load('Ba_input')
load('Sr_input')
load('Ca_input')
load('Mg_input')
load('Alk_input')

%%Defining treatment goals for contaminants (end), solubility constant
%(Ksp for Ba sulfate used), and the Ba:Sr co-precipitation ratio (barium removed
%per strontium

Ba_end = 2.62 * 10^(-5);
Sr_end = 1.80 * 10^(-2);
Ca_end = 1.22 * 10^(-1);
Mg_end = 5.76 * 10^(-2);
s = 5.67;
Ksp = 10^(-9.98);

%Setting the number of data points in the files (n) and preallocating
%vectors

n = 10000;
Sulfate = zeros(n,1);
ba_add = zeros(n,1);

%Generating a uniform distribution for the "inefficiency factor" found by
%comparing experimental data to model results

sulfate_inefficiency = uni_dist(n, 1.80, 5.77);
lime_inefficiency = uni_dist(n, 1.46, 1.88);
soda_inefficiency = uni_dist(n, 0.77, 1.80);

%Main loop for Ba and Sr removal

for i = 1:n

%Calculates change in pollutant concentrations

    Ba = Ba_input(i,1) - Ba_end;    
    if Ba < 0
        Ba = 0;
    end
    Sr = Sr_input(i,1) - Sr_end;    
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

    Sulfate(i,1) = ((Ksp / Ba_end) + Sr + Ba) .* sulfate_inefficiency(1,i);
    
%Converts lime and soda ash to g/m^3

    Sulfate(i,1) = Sulfate(i,1) * 1000 * 142.04; %Molar mass as sodium sulfate
    
end

%Main loop for Ca and Mg removal

%Setting solubility product constant for lime soda ash modeling

Ksp = 10^(-8.48);

for i = 1:n

%Calculates change in pollutant concentrations

    Ca = Ca_input(i,1) - Ca_end;
    if Ca < 0
        Ca = 0;
    end
    Mg = Mg_input(i,1) - Mg_end;    
    if Mg < 0
        Mg = 0;

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

    Lime(i,1) = CCa + 2*CMg + NCMg + (Ksp / Ca_end) .* lime_inefficiency(1,i);

    Soda(i,1) = NCCa + NCMg .* soda_inefficiency(1,i);
    
%Converts lime and soda ash to g/m^3

    Lime(i,1) = Lime(i,1) * 1000 * 100.0869;
    
    Soda(i,1) = Soda(i,1) * 1000 * 105.9888; %molar mass is for anhydrous
    
end

%Saves the sulfate, lime, and soda ash required to a file
csvwrite('Sulfate_Sulfate.csv', Sulfate);
csvwrite('Sulfate_Lime.csv', Lime);
csvwrite('Sulfate_Soda.csv', Soda);

toc