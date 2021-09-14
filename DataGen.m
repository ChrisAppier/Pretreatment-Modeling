%This script creates n data points via latin hyper cube sampling from
%empirical data

clc

%Establishing number of data points to generate
n = 30000;

%Loading alkalinity data calculated using the AlkCalc script
load('Alk_data')

%Importing the USGS water data
Ba_sample = xlsread('waterdata.xlsx','CA2:CA274');
Sr_sample = xlsread('waterdata.xlsx','ED2:ED274');
Ca_sample = xlsread('waterdata.xlsx','CH2:CH274');
Mg_sample = xlsread('waterdata.xlsx','DA2:DA274');
Na_sample = xlsread('waterdata.xlsx','DK2:DK274');
SO4_sample = xlsread('waterdata.xlsx','DV2:DV274');
Alk_sample = Alk;

%Converting the data from mg/L to mol/L
Ba_sample = Ba_sample / (1000 * 137.33);
Ca_sample = Ca_sample / (1000 * 40.078);
Sr_sample = Sr_sample / (1000 * 87.62);
Mg_sample = Mg_sample / (1000 * 24.305);
Na_sample = Na_sample / (1000 * 22.99);
SO4_sample = SO4_sample / (1000 * 96.056);

%Creating distributions based on the USGS water data
Ba_input = lhs_empir(Ba_sample, n);
Sr_input = lhs_empir(Sr_sample, n);
Ca_input = lhs_empir(Ca_sample, n);
Mg_input = lhs_empir(Mg_sample, n);
Na_input = lhs_empir(Na_sample, n);
SO4_input = lhs_empir(SO4_sample, n);
Alk_input = lhs_empir(Alk_sample, n);

%Saving the LHC distributions to matlab files
save Ba_input.mat Ba_input;
save Sr_input.mat Sr_input;
save Ca_input.mat Ca_input;
save Mg_input.mat Mg_input;
save Na_input.mat Na_input;
save SO4_input.mat SO4_input;
save Alk_input.mat Alk_input;
csvwrite('Ba_input.csv', Ba_input);
csvwrite('Sr_input.csv', Sr_input);
csvwrite('Ca_input.csv', Ca_input);
csvwrite('Mg_input.csv', Mg_input);
csvwrite('Na_input.csv', Na_input);
csvwrite('SO4_input.csv', SO4_input);
csvwrite('Alk_input.csv', Alk_input);

%Plotting the original and generated data sets
nbin = 100;

figure
subplot(2,4,1); histogram(Ba_input, nbin); title('Barium Generated'); xlabel('mg/L'); ylabel('Count');
subplot(2,4,2); histogram(Sr_input, nbin); title('Strontium Generated'); xlabel('mg/L'); ylabel('Count');
subplot(2,4,3); histogram(Ca_input, nbin); title('Calcium Generated'); xlabel('mg/L'); ylabel('Count');
subplot(2,4,4); histogram(Mg_input, nbin); title('Magnesium Generated'); xlabel('mg/L'); ylabel('Count');
subplot(2,4,5); histogram(Na_input, nbin); title('Sodium Generated'); xlabel('mg/L'); ylabel('Count');
subplot(2,4,6); histogram(SO4_input, nbin); title('Sulfate Generated'); xlabel('mg/L'); ylabel('Count');
subplot(2,4,7); histogram(Alk_input, nbin); title('Alkalinity Generated'); xlabel('mg/L'); ylabel('Count');

figure
subplot(2,4,1); histogram(Ba_sample, nbin); title('Barium Original'); xlabel('mg/L'); ylabel('Count');
subplot(2,4,2); histogram(Sr_sample, nbin); title('Strontium Original'); xlabel('mg/L'); ylabel('Count');
subplot(2,4,3); histogram(Ca_sample, nbin); title('Calcium Original'); xlabel('mg/L'); ylabel('Count');
subplot(2,4,4); histogram(Mg_sample, nbin); title('Magnesium Original'); xlabel('mg/L'); ylabel('Count');
subplot(2,4,5); histogram(Na_sample, nbin); title('Sodium Original'); xlabel('mg/L'); ylabel('Count');
subplot(2,4,6); histogram(SO4_sample, nbin); title('Sulfate Original'); xlabel('mg/L'); ylabel('Count');
subplot(2,4,7); histogram(Alk_sample, nbin); title('Alkalinity Original'); xlabel('mg/L'); ylabel('Count');





