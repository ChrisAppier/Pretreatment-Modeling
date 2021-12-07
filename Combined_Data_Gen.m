%This model takes in USGS water data, calculates alkalinity, creates a
%larger data set using LHS, calculates the limits based on RO fouling, and
%removes all water samples that require no treatment

%% Global Variables

%n = number of data points to generate, m = number of data points requiring
%treatment desired, nbin is the number of histogram bins for original vs
%generated data set plots
n = 30000;
m = 10000;
nbin = 100;

%% Alkalinity calculator

%Importing the data (mg/L) into variables
HCO3 = xlsread('waterdata.xlsx','CG2:CG274');
CO3 = xlsread('waterdata.xlsx','CF2:CF274');
pH = xlsread('waterdata.xlsx','BG2:BG274');

%Checks the number of data points
HCO3_size = size(HCO3,1);
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

%Plots the calculated alkalinity of the USGS data
plot(x, Alk, 'LineWidth', 2.0); title('Alkalinity'); ylim([0 0.1]);
xlabel('Data Point'); ylabel('Alkalinity (eq/L)');

%Stores the alkalinity
save Alk_data.mat Alk;
csvwrite('Alk_data.csv', Alk);

%% Data Generator

%Importing the USGS water data
Ba_sample = xlsread('waterdata.xlsx','CA2:CA274');
Sr_sample = xlsread('waterdata.xlsx','ED2:ED274');
Ca_sample = xlsread('waterdata.xlsx','CH2:CH274');
Mg_sample = xlsread('waterdata.xlsx','DA2:DA274');
Na_sample = xlsread('waterdata.xlsx','DK2:DK274');
SO4_sample = xlsread('waterdata.xlsx','DV2:DV274');
Alk_sample = Alk;
HCO3_sample = HCO3;

%Converting the data from mg/L to mol/L
Ba_sample = Ba_sample / (1000 * 137.33);
Ca_sample = Ca_sample / (1000 * 40.078);
Sr_sample = Sr_sample / (1000 * 87.62);
Mg_sample = Mg_sample / (1000 * 24.305);
Na_sample = Na_sample / (1000 * 22.99);
SO4_sample = SO4_sample / (1000 * 96.056);
HCO3_sample = HCO3_sample / (1000 * 61);

%Creating distributions based on the USGS water data
Ba_input = lhs_empir(Ba_sample, n);
Sr_input = lhs_empir(Sr_sample, n);
Ca_input = lhs_empir(Ca_sample, n);
Mg_input = lhs_empir(Mg_sample, n);
Na_input = lhs_empir(Na_sample, n);
SO4_input = lhs_empir(SO4_sample, n);
Alk_input = lhs_empir(Alk_sample, n);
HCO3_input = lhs_empir(HCO3_sample, n);

%Plotting the original and untrimmed generated data sets
figure
subplot(2,4,1); histogram(Ba_input, nbin); title('Barium Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,2); histogram(Sr_input, nbin); title('Strontium Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,3); histogram(Ca_input, nbin); title('Calcium Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,4); histogram(Mg_input, nbin); title('Magnesium Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,5); histogram(Na_input, nbin); title('Sodium Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,6); histogram(SO4_input, nbin); title('Sulfate Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,7); histogram(Alk_input, nbin); title('Alkalinity Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,8); histogram(HCO3_input, nbin); title('Bicarbonate Generated'); xlabel('mol/L'); ylabel('Count');

figure
subplot(2,4,1); histogram(Ba_sample, nbin); title('Barium Original'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,2); histogram(Sr_sample, nbin); title('Strontium Original'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,3); histogram(Ca_sample, nbin); title('Calcium Original'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,4); histogram(Mg_sample, nbin); title('Magnesium Original'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,5); histogram(Na_sample, nbin); title('Sodium Original'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,6); histogram(SO4_sample, nbin); title('Sulfate Original'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,7); histogram(Alk_sample, nbin); title('Alkalinity Original'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,8); histogram(HCO3_sample, nbin); title('Bicarbonate Original'); xlabel('mol/L'); ylabel('Count');

%% Limit Calculator

%SO4 Based Limits

%Converting SO4_input to mg/l
SO4_input = SO4_input * 96.056 * 1000;

%Constants
mass = SO4_input; %Mass
vol = 1; %Volume
sel = 0.95; %SO4 selectivity
rec = 0.5; %Recovery
Ksp_SO4_Ca = 6.91831*10^(-5); %Ca selectivity
Ksp_SO4_Ba = 1.0842*10^(-10); %Ba selectivity
Ksp_SO4_Sr = 3.44*10^(-7); %Sr selectivity 
Ca_ends = zeros(n,1);
Ba_ends = zeros(n,1);
Sr_ends = zeros(n,1);

%Main loop
for i=1:n
    
    perm = (1-sel)*mass(i);

    rej = (mass(i)-perm)/((1-rec)*vol); %mg/l

    rej = rej / (96*1000); %mol/l

    Ca_sat = Ksp_SO4_Ca/rej; %mol/l
    Ca_ends(i) = Ca_sat*(1-rec)/sel; %mol/l
    
    Ba_sat = Ksp_SO4_Ba/rej; %mol/l
    Ba_ends(i) = Ba_sat*(1-rec)/sel; %mol/l
        
    Sr_sat = Ksp_SO4_Sr/rej; %mol/l
    Sr_ends(i) = Sr_sat*(1-rec)/sel; %mol/l

end

%HCO3 Based Limits, only used if lower than SO4 limits

%Converting HCO3_input to mg/l
HCO3_input = HCO3_input * 1000 * 61; 

%Constants
mass = HCO3_input; %Mass
vol = 1; %Volume
sel = 0.95; %SO4 selectivity
rec = 0.5; %Recovery
Ksp_HCO3_Ca = 3.3*10^(-9); %Ca selectivity
Ksp_HCO3_Ba = 2.58*10^(-9); %Ba selectivity
Ksp_HCO3_Sr = 5.6*10^(-10); %Sr selectivity 

%Main loop
for i=1:n
    
    perm = (1-sel)*mass(i);

    rej = (mass(i)-perm)/((1-rec)*vol); %mg/l

    rej = rej / (61*1000); %mol/l

    Ca_sat = Ksp_HCO3_Ca/rej; %mol/l
    
    if Ca_sat*(1-rec)/sel < Ca_ends(i)
    Ca_ends(i) = Ca_sat*(1-rec)/sel; %mol/l
    end
    
    Ba_sat = Ksp_HCO3_Ba/rej; %mol/l
    
    if Ba_sat*(1-rec)/sel < Ba_ends(i)
    Ba_ends(i) = Ba_sat*(1-rec)/sel; %mol/l
    end
        
    Sr_sat = Ksp_HCO3_Sr/rej; %mol/l
    
    if Sr_sat*(1-rec)/sel < Sr_ends(i)
    Sr_ends(i) = Sr_sat*(1-rec)/sel; %mol/l
    end

end

%% Clean Water Removal

%l keeps track of the new data set index. Mg is not used for this because
%of calculation time in IEPSS model for generating Mg limit data
count = 0;
l = 1;

%Creating non-zero influent water vectors and their corresponding limits
Ca_input_nz = zeros(m,1);
Mg_input_nz = zeros(m,1);
Ba_input_nz = zeros(m,1);
Sr_input_nz = zeros(m,1);
Alk_input_nz = zeros(m,1);
Na_input_nz = zeros(m,1);
SO4_input_nz = zeros(m,1);
HCO3_input_nz = zeros(m,l);
Ca_ends_nz = zeros(m,1);
Mg_ends_nz = zeros(m,1);
Ba_ends_nz = zeros(m,1);
Sr_ends_nz = zeros(m,1);

%This section finds only waters that require some treatment and adds them
%into a new variable along with their corresponding limits
for i = 1:n
    if Ca_input(i) > Ca_ends(i) || Ba_input(i) > Ba_ends(i) || Sr_input(i) > Sr_ends(i)
        Ca_input_nz(l) = Ca_input(i);
        Mg_input_nz(l) = Mg_input(i);
        Ba_input_nz(l) = Ba_input(i);
        Sr_input_nz(l) = Sr_input(i);
        Na_input_nz(l) = Na_input(i);
        Alk_input_nz(l) = Alk_input(i);
        SO4_input_nz(l) = SO4_input(i);
        HCO3_input_nz(l) = HCO3_input(i);
        Ca_ends_nz(l) = Ca_ends(i);
        Ba_ends_nz(l) = Ba_ends(i);
        Sr_ends_nz(l) = Sr_ends(i);
        l = l+1;    
    end
    
    if l > m
        break;
    end
   
end

%% Renaming the generated data to match water treatment models' expected input naming system

%Resizing the original vectors to expected size/name
Ca_ends = zeros(m,1);
Ba_ends = zeros(m,1);
Sr_ends = zeros(m,1);
Ca_input = zeros(m,1);
Mg_input = zeros(m,1);
Ba_input = zeros(m,1);
Sr_input = zeros(m,1);
Alk_input = zeros(m,1);
SO4_input = zeros(m,1);
HCO3_input = zeros(m,l);

%Transferring data to original vectors
Ca_ends = Ca_ends_nz;
Ba_ends = Ba_ends_nz;
Sr_ends = Sr_ends_nz;
Ca_input = Ca_input_nz;
Mg_input = Mg_input_nz;
Ba_input = Ba_input_nz;
Sr_input = Sr_input_nz;
Na_input = Na_input_nz;
Alk_input = Alk_input_nz;
SO4_input = SO4_input_nz;
HCO3_input = HCO3_input_nz;

%Saving the data (mol/l)
save Ba_input.mat Ba_input;
save Sr_input.mat Sr_input;
save Ca_input.mat Ca_input;
save Mg_input.mat Mg_input;
save Na_input.mat Na_input;
save SO4_input.mat SO4_input;
save Alk_input.mat Alk_input;
save HCO3_input.mat HCO3_input;
save Ca_ends.mat Ca_ends;
save Ba_ends.mat Ba_ends;
save Sr_ends.mat Sr_ends;
csvwrite('Ba_input.csv', Ba_input);
csvwrite('Sr_input.csv', Sr_input);
csvwrite('Ca_input.csv', Ca_input);
csvwrite('Mg_input.csv', Mg_input);
csvwrite('Na_input.csv', Na_input);
csvwrite('SO4_input.csv', SO4_input);
csvwrite('Alk_input.csv', Alk_input);

%Plots the trimmed generated data
figure
subplot(2,4,1); histogram(Ba_input, nbin); title('Trimmed Barium Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,2); histogram(Sr_input, nbin); title('Trimmed Strontium Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,3); histogram(Ca_input, nbin); title('Trimmed Calcium Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,4); histogram(Mg_input, nbin); title('Trimmed Magnesium Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,5); histogram(Na_input, nbin); title('Trimmed Sodium Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,6); histogram(SO4_input, nbin); title('Trimmed Sulfate Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,7); histogram(Alk_input, nbin); title('Trimmed Alkalinity Generated'); xlabel('mol/L'); ylabel('Count');
subplot(2,4,8); histogram(HCO3_input, nbin); title('Trimmed Bicarbonate Generated'); xlabel('mol/L'); ylabel('Count');

%Debug
disp('Ca Limit (mol/L): ') 
disp(mean(Ca_ends))
disp('Ba Limit (mol/L): ') 
disp(mean(Ba_ends))
disp('Sr Limit (mol/L): ') 
disp(mean(Sr_ends))
