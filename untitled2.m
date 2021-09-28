HCO3 = xlsread('waterdata.xlsx','CG2:CG274');
CO3 = xlsread('waterdata.xlsx','CF2:CF274');

disp(mean(HCO3))
disp(mean(CO3))