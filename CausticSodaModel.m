%This model calculates the amount of caustic soda required to remove
%calcium, magnesium, barium, and strontium hardness in water.

clc
clear 
tic

%Loading the input data files as mol/L

load('Ca_input')
load('Mg_input')
load('Ba_input')
load('Sr_input')
load('Alk_input')
load('inefficiency')

%Defining treatment goals for contaminants (mol/L), solubility constant
%(Ksp for Ca carbonate used), and the kWh energy required per mol/L of caustic soda (e_fac).

Ca_end = 4.37 * 10^(-4);
Mg_end = 8.23 * 10^(-4);
Ba_end = 4.00 * 10^(-4);
Sr_end = 3.14 * 10^(-4);
Ksp = 10^(-8.48);
e_fac  = 0.159988;
    
%Setting the number of data points in the files (n) and preallocating
%vectors

n = 5;
Caustic = zeros(n,1);
Soda = zeros(n,1);
Energy = zeros(n,1);

%Generating a uniform distribution for the "inefficiency factor" found by
%comparing experimental data to model results

caustic_inefficiency = uni_dist(n, 1.64, 2.13);
soda_inefficiency = uni_dist(n, 0.77, 1.80);

%Main loop

for i = 1:n

%Calculates change in pollutant concentrations

    Ca = Ca_input(1,i) - Ca_end;

    Mg = Mg_input(1,i) - Mg_end;

    Ba = Ba_input(1,i) - Ba_end;

    Sr = Sr_input(1,i) - Sr_end;

%Calculates carbonate hardness (CCa and CMg) and non-carbonate hardness
%(NCCa and NCMg). Ba and Sr behave similar to Ca and are all combined under
%the CCa and NCCa variables

    if 2 * (Ca + Ba + Sr + Mg) <= Alk_input(1,i)
    
        CCa = Ca + Ba + Sr;
        NCCa = 0;
        CMg = Mg;
        NCMg = 0;

    elseif 2 * (Ca + Ba + Sr + Mg) > Alk_input(1,i) && 2 * (Ca + Ba + Sr) <= Alk_input(1,i)
    
        CCa = Ca + Ba + Sr;
        NCCa = 0;
        CMg = (Alk_input(1,i) - 2 * (Ca + Ba + Sr)) / 2;
        NCMg = Mg - CMg;
    
    elseif 2* (Ca + Ba + Sr + Mg) > Alk_input(1,i) && 2* (Ca + Ba + Sr) > Alk_input(1,i)
    
        CCa = Alk_input(1,i);
        NCCa = Ca + Ba + Sr - (0.5 * Alk_input(1,i));
        CMg = 0;
        NCMg = Mg;
    
    end

%Calculates the amount of caustic soda, soda ash required, and energy
%required

    Caustic(i,1) = 2*CCa + 4*CMg + 2*NCMg + (Ksp / Ca_end) .* caustic_inefficiency(1,i);

    Soda(i,1) = NCCa .* soda_inefficiency(1,i);
    
    Energy(i,1) = Caustic(i,1) * e_fac;
    
%Converts caustic soda and soda ash to g/m^3

    Caustic(i,1) = Caustic(i,1) * 1000 * 39.9971;
    
    Soda(i,1) = Soda(i,1) * 1000 * 105.9888; %molar mass is for anhydrous

end
    
%Saves the energy, caustic soda, and soda ash required to a file
csvwrite('CausticSoda_Caustic.csv', Caustic);
csvwrite('CausticSoda_Soda.csv', Soda);
csvwrite('CausticSoda_Energy.csv', Energy);

toc