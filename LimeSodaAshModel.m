%This model calculates the amount of lime and soda ash required to remove
%calcium, magnesium, barium, and strontium hardness in water. Strontium
%removal is through substitution with calcium.

%Loading the input data files as mol/L

load('Ca_input')
load('Mg_input')
load('Ba_input')
load('Sr_input')
load('Alk_input')
load('inefficiency')

%Defining treatment goals for contaminants (end), solubility constant
%(Ksp for Ca used), and the Ca:Sr co-precipitation ratio (calcium removed
%per strontium)

Ca_end = 0.01;
Mg_end = 0.01;
Ba_end = 0.01;
Sr_end = 0.01;
s = 3;
Ksp = 10 * 10^(-8);

%Checking the number of data points for each input to ensure they are equal and
%setting the number of loops for the model equal to that size or ending the
%program. Also preallocates lime and soda ash addition variables for
%efficiency.

if isequal(height(Ca_input),height(Mg_input),height(Ba_input),height(Sr_input),height(inefficiency),height(Alk_input))
    
    n = height(Ba_input);
    
    Lime = zeros(n,1);
    
    Soda = zeros(n,1);
    
else

    disp('Input files contain unequal number of data points. The program will not run. Check the data and try again')      
    return
    
end

%Main loop

for i = 1:n

%Calculates change in pollutant concentrations

    Ca = Ca_input(i) - Ca_end;

    Mg = Mg_input(i) - Mg_end;

    Ba = Ba_input(i) - Ba_end;

    Sr = Sr_input(i) - Sr_end;

%Checking if enough Ca is present to remove Sr and adds soda ash if not.
%First if statement is for all calcium carbonate hardness, second if
%statement is for combination of carbonate and non-carbonate calcium
%hardness

    if 2 * Ca < s * Sr && 2 * Ca <= Alk_input(i)
    
        Ca_add = (s * Sr) - (2 * Ca);
    
        Ca = Ca + Ca_add;
    
    end

    if 2 * Ca < s * Sr && 2 * Ca > Alk_input(i)
    
        Ca_add = (s * Sr) - Alk_input(i) - (Ca - (0.5 * Alk_input(i));
    
        Ca = Ca + Ca_add;
    
    end

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
    
    elseif 2 * (Ca + Ba + Sr + Mg) > Alk_input(i) && 2 * (Ca + Ba + Sr) > Alk_input(i)
    
        CCa = Alk_input(i);
        NCCa = Ca + Ba + Sr - (0.5 * Alk_input(i));
        CMg = 0;
        NCMg = Mg;
    
    end

%Calculates the amount of lime and soda ash required

    Lime(i) = CCa + 2*CMg + NCMg + (Ksp / Ca_end);

    Soda(i) = NCCa + NCMg + Ca_add;
    
end

%Transposes the columns to rows and writes the lime and soda ash additions to a text file
Lime = transpose(Lime);
Soda = transpose(Soda);
writematrix(Lime,'Lime_Addition.txt');
writematrix(Soda, 'Lime_SodaAsh_Addition.txt');