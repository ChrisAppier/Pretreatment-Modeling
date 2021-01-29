%This script creates n data points via latin hyper cube sampling from
%empirical data

clc

%Importing the USGS water data
Ba_sample = xlsread('waterdata.xlsx','E3:E447');
Sr_sample = xlsread('waterdata.xlsx','BI3:BI488');
Ca_sample = xlsread('waterdata.xlsx','K3:K2758');
Mg_sample = xlsread('waterdata.xlsx','AQ3:AQ2740');
Na_sample = xlsread('waterdata.xlsx','AU3:AU2614');
Alk_sample = xlsread('waterdata.xlsx','A3:A7');

%Creating distributions based on the USGS water data
Ba_input = lhs_empir(Ba_sample, 10000);
Sr_input = lhs_empir(Sr_sample, 10000);
Ca_input = lhs_empir(Ca_sample, 10000);
Mg_input = lhs_empir(Mg_sample, 10000);
Na_input = lhs_empir(Na_sample, 10000);
Alk_input = lhs_empir(Alk_sample, 10000);

%Saving the LHC distributions to matlab files
save Ba_input.mat Ba_input;
save Sr_input.mat Sr_input;
save Ca_input.mat Ca_input;
save Mg_input.mat Mg_input;
save Na_input.mat Na_input;
save Alk_input.mat Alk_input;


%Plotting the original and generated data sets
nbin = 100;

figure
subplot(2,3,1); histogram(Ba_input, nbin); title('Barium Generated');
subplot(2,3,2); histogram(Sr_input, nbin); title('Strontium Generated');
subplot(2,3,3); histogram(Ca_input, nbin); title('Calcium Generated');
subplot(2,3,4); histogram(Mg_input, nbin); title('Magnesium Generated');
subplot(2,3,5); histogram(Na_input, nbin); title('Sodium Generated');
subplot(2,3,6); histogram(Alk_input, nbin); title('Alkalinity Generated');

figure
subplot(2,3,1); histogram(Ba_sample, nbin); title('Barium Original');
subplot(2,3,2); histogram(Sr_sample, nbin); title('Strontium Original');
subplot(2,3,3); histogram(Ca_sample, nbin); title('Calcium Original');
subplot(2,3,4); histogram(Mg_sample, nbin); title('Magnesium Original');
subplot(2,3,5); histogram(Na_sample, nbin); title('Sodium Original');
subplot(2,3,6); histogram(Alk_sample, nbin); title('Alkalinity Original');