%This model calculates the amount of caustic soda required to remove
%calcium, magnesium, barium, and strontium hardness in water.

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

%Uses the solubility constant that is the highest (Ba)

Ksp = 5.0 * 10^(-3);

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

Caustic = 2*CCa + 4*CMg + 2*NCMg + (Ksp / Ba_end);

Soda = NCCa;
    
%Debug Outputs
%disp('CCa =');
%disp(CCa);
%disp('NCCa =');
%disp(NCCa);
%disp('CMg');
%disp(CMg);
%disp('NCMg =');
%disp(NCMg);

%Additions required output
disp('Caustic Soda =');
disp(Caustic);
disp('Soda Ash =')
disp(Soda);