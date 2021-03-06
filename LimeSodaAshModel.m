%This model calculates the amount of lime and soda ash required to remove
%calcium, magnesium, barium, and strontium hardness in water. Strontium
%removal is through substitution with calcium.

clc
clear 
tic

%Loading the input data files as mol/L

load('Ca_input')
load('Mg_input')
load('Ba_input')
load('Sr_input')
load('Alk_input')

%Defining treatment goals for contaminants (end), solubility constant
%(Ksp for Ca carbonate used)

Ba_end = 0 * 10^(-5);
Ca_end = 6.363 * 10^(-4);
Mg_end = 8.091 * 10^(-4);
Ksp = 10^(-8.48);

%Setting the number of data points in the files (n) and preallocating
%vectors

n = 10000;
Lime = zeros(n:1);
Soda = zeros(n:1);
Sr_end = zeros(n:1);

%Generating a uniform distribution for the "inefficiency factor" found by
%comparing experimental data to model results

lime_inefficiency = uni_dist(n, 1.46, 1.88);
soda_inefficiency = uni_dist(n, 0.77, 1.80);

%Main loop

for i = 1:n

%Calculates change in pollutant concentrations

    Ca = Ca_input(i,1) - Ca_end;
    if Ca < 0
        Ca = 0;
    end
    Mg = Mg_input(i,1) - Mg_end;    
    if Mg < 0
        Mg = 0;
    end
    Ba = Ba_input(i,1) - Ba_end;    
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

    Lime(i,1) = CCa + 2*CMg + NCMg + (Ksp / Ca_end) .* lime_inefficiency(1,i);

    Soda(i,1) = NCCa + NCMg .* soda_inefficiency(1,i);
    
%Converts lime and soda ash to g/m^3

    Lime(i,1) = Lime(i,1) * 1000 * 100.0869;
    
    Soda(i,1) = Soda(i,1) * 1000 * 105.9888; %molar mass is for anhydrous
    
%Calculates the amount of Sr left in water
    
    Sr_end(i,1) = Sr_input(i,1) - Sr;
    
end

%Saves the lime and soda ash required to a file
csvwrite('LimeSodaAsh_Lime.csv', Lime);
csvwrite('LimeSodaAsh_Soda.csv', Soda);

toc