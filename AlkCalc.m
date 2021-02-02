%This script takes water data from an Excel spreadsheet, calculates the
%alkalinity, and inserts the alkalinity data into a stored variable for
%other scripts to use

%Importing the data into variables
HCO3 = xlsread('waterdata.xlsx','A1:A1');
CO3 = xlsread('waterdata.xlsx','B1:B1');
pH = xlsread('waterdata.xlsx','C1:C1');

%Checks the number of data points
HCO3_size = size(HCO3, 1); %may need to swap 1 and vector%
CO3_size = size(CO3, 1);
pH_size = size(pH, 1);

%Ends the script if the sizes of the variables differ
if HCO3_size ~= CO3_size || HCO3_size ~= pH_size
    disp('Error: Variable sizes are not the same')
    return
end

%

%Calculates the alkalinity
parfor i = 1:pH_size
    
   Alk = ; %find equation to fill in the rest 
    
end

%Stores the alkalinity
save Alk_data.mat Alk;