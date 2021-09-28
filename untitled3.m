SO4_sample = xlsread('waterdata.xlsx','DV2:DV274');
Ca_sample = xlsread('waterdata.xlsx','CH2:CH274');

disp(mean(Ca_sample))
disp(median(SO4_sample))