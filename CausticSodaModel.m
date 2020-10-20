%This model calculates the amount of caustic soda required to remove
%calcium, magnesium, barium, and strontium hardness in water.

%Loading the input data files as mol/L

load('Ca_input')
load('Mg_input')
load('Ba_input')
load('Sr_input')
load('Alk_input')
load('inefficiency')

%Defining constants, these still need to be filled in with real #s

Ca_end = 0.01;
Mg_end = 0.01;
Ba_end = 0.01;
Sr_end = 0.01;
Ksp = 10 * 10^(-8);

%Checking the number of data points for each input to ensure they are equal and
%setting the number of loops for the model equal to that size or ending the
%program

if isequal(height(Ca_input),height(Mg_input),height(Ba_input),height(Sr_input),height(inefficiency),height(Alk_input))
    
    n = height(Ba_input);
    
else

    disp('Input files contain unequal number of data points. The program will not run. Check the data and try again')      
    return
    
end

%Main loop

for i = 1:n

%Calculates change in pollutant concentrations

    Ca = Ca_start - Ca_end;

    Mg = Mg_start - Mg_end;

    Ba = Ba_start - Ba_end;

    Sr = Sr_start - Sr_end;

%Calculates carbonate hardness (CCa and CMg) and non-carbonate hardness
%(NCCa and NCMg). Ba and Sr behave similar to Ca and are all combined under
%the CCa and NCCa variables

    if 2 * (Ca + Ba + Sr + Mg) <= Alk_input(i)
    
        CCa = Ca + Ba + Sr;
        NCCa = 0;
        CMg = Mg;
        NCMg = 0;

    elseif 2 * (Ca + Ba + Sr + Mg) > Alk_input(i) && 2 * (Ca + Ba + Sr) <= Alk_input(i)
    
        CCa = Ca + Ba + Sr;
        NCCa = 0;
        CMg = (Alk_input(i) - 2 * (Ca + Ba + Sr)) / 2;
        NCMg = Mg - CMg;
    
    elseif 2* (Ca + Ba + Sr + Mg) > Alk_input(i) && 2* (Ca + Ba + Sr) > Alk_input(i)
    
        CCa = Alk_input(i);
        NCCa = Ca + Ba + Sr - (0.5 * Alk_input(i));
        CMg = 0;
        NCMg = Mg;
    
    end

%Calculates the amount of caustic soda and soda ash required

    Caustic(i) = 2*CCa + 4*CMg + 2*NCMg + (Ksp / Ba_end);

    Soda(i) = NCCa;

end
    
%Transposes the columns to rows and writes the caustic soda and soda ash additions to a text file
Caustic = transpose(Caustic);
Soda = transpose(Soda);
writematrix(Caustic,'CausticSoda_Addition.txt');
writematrix(Soda, 'SodaAsh_Addition.txt');