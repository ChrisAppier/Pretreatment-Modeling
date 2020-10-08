%This model calculates the amount of sodium sulfate required to remove
%calcium, barium, and strontium in water. Strontium removal is through
%substitution with barium.

%Prompting the user for the starting and ending Ca, Ba, and Sr concentrations

prompt = 'Enter calcium concentration (mol/L)\n';
Ca_start = input(prompt);

prompt = 'Enter barium concentration (mol/L)\n';
Ba_start = input(prompt);

prompt = 'Enter strontium concentration (mol/L)\n';
Sr_start = input(prompt);

prompt = 'Enter desired calcium concentration (mol/L)\n';
Ca_end = input(prompt);

prompt = 'Enter desired barium concentration (mol/L)\n';
Ba_end = input(prompt);

prompt = 'Enter desired strontium concentration (mol/L)\n';
Sr_end = input(prompt);

%Establishes substitution ratio (Ba per Sr substitutions) and BaSO4 Ksp

s = 3;

Ksp = 9.1 * 10^(-6);

%Calculates change in pollutant concentrations

Ca = Ca_start - Ca_end;

Ba = Ba_start - Ba_end;

Sr = Sr_start - Sr_end;

%Determines if Ba will limit Sr removal

if Ba < s * (Sr)
    
    Ba_add = (s * Sr) - Ba;
    
    Ba = Ba + Ba_add;
    
end    

%Calculates the sulfate required to remove the Ba/Sr and the amount of Ba
%addition required, if necessary

sulfate = (Ksp / Ca_end) + Sr + Ba + Ca;
    
%Debug Section
%disp('Ba');
%disp(Ba);

%Additions required output
disp('Sulfate (mol/L) =');
disp(sulfate);
disp('Ba to add (mol/L) =')
disp(Ba_add);