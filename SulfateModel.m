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
load('inefficiency')

%%Defining treatment goals for contaminants (end), solubility constant
%(Ksp for Ba sulfate used), and the Ba:Sr co-precipitation ratio (barium removed
%per strontium

Ba_end = 2.62 * 10^(-5);
Sr_end = 1.80 * 10^(-2);
Ca_end = 0.00 * 10^(-5);
Mg_end = 3.20 * 10^(-4);
s = 5.67;
Ksp = 10^(-9.98);

%Setting the number of data points in the files (n) and preallocating
%vectors

n = 2;
Sulfate = zeros(n,1);
ba_add = zeros(n,1);

%Main loop for Ba and Sr removal

parfor i = 1:n

%Calculates change in pollutant concentrations

    Ba = Ba_input(1,i) - Ba_end;

    Sr = Sr_input(1,i) - Sr_end;

%Determines if Ba will limit Sr removal

    if Ba < s * Sr
    
        Ba_add(i,1) = (s * Sr) - Ba;
    
        Ba = Ba + Ba_add(i,1);
      
    else
        
        Ba_add(i,1) = 0;
    
    end    

%Calculates the sulfate required to remove the Ba and Sr

    Sulfate(i,1) = ((Ksp / Ba_end) + Sr + Ba) .* inefficiency(1,i);
    
%Converts lime and soda ash to g/m^3

    Sulfate(i,1) = Sulfate(i,1) * 1000 * 142.04; %Molar mass as sodium sulfate
    
end

%Main loop for Ca and Mg removal

%Setting solubility product constant for lime soda ash modeling

Ksp = 10^(-8.48);

parfor i = 1:n

%Calculates change in pollutant concentrations

    Ca = Ca_input(1,i) - Ca_end;

    Mg = Mg_input(1,i) - Mg_end;

%Calculates carbonate hardness (CCa and CMg) and non-carbonate hardness
%(NCCa and NCMg).

    if 2 * (Ca + Mg) <= Alk_input(1,i)
    
        CCa = Ca;
        NCCa = 0;
        CMg = Mg;
        NCMg = 0;

    elseif 2 * (Ca + Mg) > Alk_input(1,i) && 2 * Ca <= Alk_input(1,i)
    
        CCa = Ca;
        NCCa = 0;
        CMg = (Alk_input(1,i) - 2 * Ca) / 2;
        NCMg = Mg - CMg;
    
    elseif 2 * (Ca + Mg) > Alk_input(1,i) && 2 * Ca > Alk_input(1,i)
    
        CCa = Alk_input(1,i);
        NCCa = Ca - (0.5 * Alk_input(1,i));
        CMg = 0;
        NCMg = Mg;
    
    end

%Calculates the amount of lime and soda ash required

    Lime(i,1) = CCa + 2*CMg + NCMg + (Ksp / Ca_end);

    Soda(i,1) = NCCa + NCMg;
    
%Converts lime and soda ash to g/m^3

    Lime(i,1) = Lime(i,1) * 1000 * 100.0869;
    
    Soda(i,1) = Soda(i,1) * 1000 * 105.9888; %molar mass is for anhydrous
    
end

%Add lime soda ash removal model for Ca, Mg

%Saves the sulfate, lime, and soda ash required to a file
csvwrite('Sulfate_Sulfate.csv', Sulfate);
csvwrite('Sulfate_Lime.csv', Lime);
csvwrite('Sulfate_Soda.csv', Soda);

toc