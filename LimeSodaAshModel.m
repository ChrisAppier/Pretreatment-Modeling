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

n = 2;
Lime = zeros(n:1);
Soda = zeros(n:1);
Sr_end = zeros(n:1);

%Generating a uniform distribution for the "inefficiency factor" found by
%comparing experimental data to model results

lime_inefficiency = uni_dist(n, 1.46, 1.88);
soda_inefficiency = uni_dist(n, ); %Need to find #s for soda ash!!

%Main loop

for i = 1:n

%Calculates change in pollutant concentrations

    Ca = Ca_input(1,i) - Ca_end;

    Mg = Mg_input(1,i) - Mg_end;

    Ba = Ba_input(1,i) - Ba_end;

    Sr = Sr_input(1,i) * (Ca / Ca_input(1,i));


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
    
    elseif 2 * (Ca + Ba + Sr + Mg) > Alk_input(1,i) && 2 * (Ca + Ba + Sr) > Alk_input(1,i)
    
        CCa = Alk_input(1,i);
        NCCa = Ca + Ba + Sr - (0.5 * Alk_input(1,i));
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
    
    Sr_end(i,1) = Sr_input(1,i) - Sr;
    
end

%Saves the lime and soda ash required to a file
csvwrite('LimeSodaAsh_Lime.csv', Lime);
csvwrite('LimeSodaAsh_Soda.csv', Soda);

toc