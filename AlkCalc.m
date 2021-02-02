%This script takes water data from an Excel spreadsheet, calculates the
%alkalinity, and inserts the alkalinity data into a stored variable for
%other scripts to use.

%Importing the data (mg/L) into variables
HCO3 = xlsread('waterdata.xlsx','A2:A4');
CO3 = xlsread('waterdata.xlsx','B2:B4');
pH = xlsread('waterdata.xlsx','C2:C4');

%Checks the number of data points
HCO3_size = size(HCO3,1); %may need to swap 1 and vector%
CO3_size = size(CO3,1);
pH_size = size(pH,1);

%Ends the script if the sizes of the variables differ
if HCO3_size ~= CO3_size || HCO3_size ~= pH_size
    disp('Error: Variable sizes are not the same')
    return
end

%Converting the data to mol/L
HCO3 = HCO3 / (61.016 * 1000);
CO3 = CO3 / (60.008 * 1000);

%Predefining variables
pOH = zeros(size(pH,1),1);

%Calculates the H and OH concentrations based on pH
parfor i = 1:pH_size
   
    pOH = 14 - pH(i);
    
    H(i) = 10^(-pH(i));
    
    OH(i) = 10^(-pOH(i));
    
end

%Calculates the alkalinity
parfor i = 1:pH_size
    
   Alk = HCO3(i) + 2*CO3(i) + OH(i) - H(i);
    
end

%Stores the alkalinity
save Alk_data.mat Alk;