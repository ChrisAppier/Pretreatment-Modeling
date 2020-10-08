%This model calculates the amount of lime and soda ash required to remove
%calcium, magnesium, barium, and strontium hardness in water. Strontium
%removal is through substitution with calcium.

prompt = 'Enter calcium concentration (mol/L)\n';
Ca_start = input(prompt);

prompt = 'Enter magnesium concentration (mol/L)\n';
Mg_start = input(prompt);

prompt = 'Enter barium concentration (mol/L)\n';
Ba_start = input(prompt);

prompt = 'Enter strontium concentration (mol/L)\n';
Sr_start = input(prompt);

prompt = 'Enter the alkalinity (mol/L)\n';
Alk = input(prompt);

prompt = 'Enter desired calcium concentration (mol/L)\n';
Ca_end = input(prompt);

prompt = 'Enter desired magnesium concentration (mol/L)\n';
Mg_end = input(prompt);

prompt = 'Enter desired barium concentration (mol/L)\n';
Ba_end = input(prompt);

prompt = 'Enter desired strontium concentration (mol/L)\n';
Sr_end = input(prompt);

%Calculates change in pollutant concentrations

Ca = Ca_start - Ca_end;

Mg = Mg_start - Mg_end;

Ba = Ba_start - Ba_end;

Sr = Sr_start - Sr_end;

%Declaires Sr to Ca removal ratio (Ca removed per Sr) and CaCO3 Ksp

s = 3;

Ksp = 3.36 * 10^(-9);

%Checking if enough Ca is present to remove Sr and adds soda ash if not
%(needs fixed to account for NCCa)

if 2 * Ca < s * Sr
    
    Ca_add = (s * Sr) - (2 * Ca);
    
    Ca = Ca + Ca_add;
    
end

%Calculates carbonate hardness (CCa and CMg) and non-carbonate hardness
%(NCCa and NCMg). Ba and Sr behave similar to Ca and are all combined under
%the CCa and NCCa variables

if 2 * (Ca + Ba + Sr + Mg) <= Alk
    
    CCa = Ca + Ba + Sr;
    NCCa = 0;
    CMg = Mg;
    NCMg = 0;

elseif 2 * (Ca + Ba + Sr + Mg) > Alk && 2 * (Ca + Ba + Sr) <= Alk
    
    CCa = Ca + Ba + Sr;
    NCCa = 0;
    CMg = (Alk - 2 * (Ca + Ba + Sr)) / 2;
    NCMg = Mg - CMg;
    
elseif 2* (Ca + Ba + Sr + Mg) > Alk && 2* (Ca + Ba + Sr) > Alk
    
    CCa = Alk;
    NCCa = Ca + Ba + Sr - (0.5 * Alk);
    CMg = 0;
    NCMg = Mg;
    
end

%Calculates the amount of lime and soda ash required

Lime = CCa + 2*CMg + NCMg + (Ksp / Ca_end);

Soda = NCCa + NCMg + Ca_add;

%Debug Section
%disp('CCa =');
%disp(CCa);
%disp('NCCa =');
%disp(NCCa);
%disp('CMg');
%disp(CMg);
%disp('NCMg =');
%disp(NCMg);

%Additions required output
disp('Lime =');
disp(Lime);
disp('Soda Ash =')
disp(Soda);