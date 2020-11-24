%This model calculates the amount of sodium sulfate required to remove
%barium and strontium in water. Strontium removal is through
%substitution with barium.

clc 
clear

%Loading the input data files as mol/L

load('Ba_input')
load('Sr_input')
load('inefficiency')

%%Defining treatment goals for contaminants (end), solubility constant
%(Ksp for Ba sulfate used), and the Ba:Sr co-precipitation ratio (barium removed
%per strontium

Ba_end = 1.46 * 10^(-5);
Sr_end = 1.71 * 10^(-5);
s = 5.67;
Ksp = 1.1 * 10^(-10);

%Setting the number of data points in the files (n) and preallocating
%vectors

n = 3;
sulfate = zeros(n,1);
ba_add = zeros(n,1);

%Main loop

for i = 1:n

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

    sulfate(i,1) = ((Ksp / Ba_end) + Sr + Ba) .* inefficiency(1,i);
    
end

%Add lime soda ash removal model for Ca, Mg

%Transposes the columns to rows and writes the sulfate and barium additions to a text file
writematrix(sulfate,'Sulfate_Addition.txt');
writematrix(Ba_add, 'Sulfate_Ba_Addition.txt');

%Convert to mass (to match SimaPro output)