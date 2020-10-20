%This model calculates the amount of sodium sulfate required to remove
%barium and strontium in water. Strontium removal is through
%substitution with barium.

%Loading the input data files as mol/L

load('Ba_input')
load('Sr_input')
load('inefficiency')

%Defining constants, these still need to be filled in with real #s

Ba_end = 0.01;
Sr_end = 0.01;
s = 3;
Ksp = 10 * 10^(-8);

%Checking the number of data points for each input to ensure they are equal and
%setting the number of loops for the model equal to that size or ending the
%program

if isequal(height(Ba_input),height(Sr_input),height(inefficiency))
    
    n = height(Ba_input);
    
else

    disp('Input files contain unequal number of data points. The program will not run. Check the data and try again')      
    return
    
end

%Main loop

for i = 1:n

%Calculates change in pollutant concentrations

    Ba = Ba_input(i) - Ba_end;

    Sr = Sr_input(i) - Sr_end;

%Determines if Ba will limit Sr removal

    if Ba < s * Sr
    
        Ba_add(i) = (s * Sr) - Ba;
    
        Ba = Ba + Ba_add(i);
      
    else
        
        Ba_add(i) = 0;
    
    end    

%Calculates the sulfate required to remove the Ba and Sr

    sulfate(i) = ((Ksp / Ba_end) + Sr + Ba) * inefficiency(i);
    
end
    
%Debug Section


%Additions required output
%disp('Sulfate (mol/L) =');
%disp(sulfate);
%disp('Ba to add (mol/L) =')
%disp(Ba_add);