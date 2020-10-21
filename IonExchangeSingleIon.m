%This model calculates the amount of Ba removed through ion exchange resin

%Loading the input data files as mol/L

load('Ba_input')
load('Na_input')
load('inefficiency')

%Defining constants, these still need to be filled in with real #s

Ba_end = 0.01;
K = 10;
regen_percent = 90;

%Checking the number of data points for each input to ensure they are equal and
%setting the number of loops for the model equal to that size or ending the
%program. Also preallocates resin and salt variables for efficiency.

if isequal(height(Ba_input),height(Na_input),height(inefficiency))
    
    n = height(Ba_input);
    
    resin = zeros(n,1);
    
    salt = zeros(n,1);
    
else

    disp('Input files contain unequal number of data points. The program will not run. Check the data and try again')      
    return
    
end

%Main loop

for i = 1:n
    
%Calculates intermediate values to simplify resin calculations

    RBa = Ba_input(i) - Ba_end;

    Na = Na_input(i) + RBa;
    
    RNa = sqrt(((Na^2 + RBa)^2) * RBa / (K * Ba_end));

%Calculates the amount of resin required

    resin(i) = RNa + (2 * (K * (RNa^2) * Ba_end) / (Na^2));

%Calculates the salt concentration required for regeneration of the resin
%to a certain percent

    salt(i) = 0;
    
end
    
%Transposes the columns to rows and writes the sulfate and barium additions to a text file
resin = transpose(resin);
salt = transpose(salt);
writematrix(resin,'Resin_Addition.txt');
writematrix(salt, 'Salt_Addition.txt');