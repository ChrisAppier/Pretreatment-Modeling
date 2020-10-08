%This model calculates the amount of lime and soda ash required to remove
%calcium and magnesium hardness in water

prompt = 'Enter calcium concentration (mol/L)\n';
Ca_start = input(prompt);

prompt = 'Enter magnesium concentration (mol/L)\n';
Mg_start = input(prompt);

prompt = 'Enter the alkalinity (mol/L)\n';
Alk = input(prompt);

prompt = 'Enter desired calcium concentration (mol/L)\n';
Ca_end = input(prompt);

prompt = 'Enter desired magnesium concentration (mol/L)\n';
Mg_end = input(prompt);

%Calculates change in pollutant concentrations

Ca = Ca_start - Ca_end;

Mg = Mg_start - Mg_end;

%Calculates carbonate hardness (CCa and CMg) and non-carbonate hardness
%(NCCa and NCMg)

if 2 * (Ca + Mg) <= Alk
    
    CCa = Ca;
    NCCa = 0;
    CMg = Mg;
    NCMg = 0;

elseif 2 * (Ca + Mg) > Alk && 2 * Ca <= Alk
    
    CCa = Ca;
    NCCa = 0;
    CMg = Alk - (2 * Ca);
    NCMg = Mg - CMg;
    
elseif 2 * (Ca + Mg) > Alk && 2 * Ca > Alk
    
    CCa = Alk / 2;
    NCCa = Ca - CCa;
    CMg = 0;
    NCMg = Mg;
    
end    

%Calculates the amount of lime and soda ash required

Lime = CCa + 2*CMg + NCMg + (Ksp / Ca_end);

Soda = NCCa + NCMg;
    
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
disp('Lime =');
disp(Lime);
disp('Soda Ash =')
disp(Soda);