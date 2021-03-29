%This script takes water data from an Excel spreadsheet in mg/L, calculates the
%alkalinity in eq/L, and inserts the alkalinity data into a stored variable for
%other scripts to use.

%Importing the data (mg/L) into variables
HCO3 = xlsread('waterdata.xlsx','CG2:CG274');
CO3 = xlsread('waterdata.xlsx','CF2:CF274');
pH = xlsread('waterdata.xlsx','BG2:BG274');

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
Alk = zeros(size(pH,1),1);
x = zeros(size(pH,1),1);

%Calculates the H and OH concentrations based on pH
for i = 1:pH_size
   
    pOH(i) = 14 - pH(i);
    
    H(i) = 10^(-pH(i));
    
    OH(i) = 10^(-pOH(i));
    
end

%Calculates the alkalinity
for i = 1:pH_size
    
   Alk(i) = HCO3(i) + 2*CO3(i) + OH(i) - H(i);
    
end

for i = 1:pH_size
    x(i) = i;
end

plot(x, Alk, 'LineWidth', 2.0); title('Alkalinity'); ylim([0 0.1]);
xlabel('Data Point'); ylabel('Alkalinity (eq/L)');

%Stores the alkalinity
save Alk_data.mat Alk;
csvwrite('Alk_data.csv', Alk);