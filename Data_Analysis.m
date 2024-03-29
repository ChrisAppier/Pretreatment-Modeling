%This script reads the data generated by the treatment models, analyzes it,
%and outputs some statistics and graphs

%Loading input/output values
load('Ca_input')
load('Mg_input')
load('Ba_input')
load('Sr_input')
load('Ca_ends')
%load('Mg_ends')
load('Ba_ends')
load('Sr_ends')

%Loading IEPSS model outputs
IEPSS_AP = csvread('IEPSS_AP.csv');
IEPSS_EP = csvread('IEPSS_EP.csv');
IEPSS_GWP = csvread('IEPSS_GWP.csv');
IEPSS_ODP = csvread('IEPSS_ODP.csv');
IEPSS_POCP = csvread('IEPSS_POCP.csv');
IEPSS_PEU = csvread('IEPSS_PEU.csv');
IEPSS_CAR = csvread('IEPSS_CAR.csv');
IEPSS_NCAR = csvread('IEPSS_NCAR.csv');
IEPSS_RES = csvread('IEPSS_RES.csv');
IEPSS_ETX = csvread('IEPSS_ETX.csv');
%IEPSS_Mg_ends = csvread('Mg_ends.csv');

%Loading lime soda ash model outputs
LimeSodaAsh_AP = csvread('LimeSodaAsh_AP.csv');
LimeSodaAsh_EP = csvread('LimeSodaAsh_EP.csv');
LimeSodaAsh_GWP = csvread('LimeSodaAsh_GWP.csv');
LimeSodaAsh_ODP = csvread('LimeSodaAsh_ODP.csv');
LimeSodaAsh_POCP = csvread('LimeSodaAsh_POCP.csv');
LimeSodaAsh_PEU = csvread('LimeSodaAsh_PEU.csv');
LimeSodaAsh_CAR = csvread('LimeSodaAsh_CAR.csv');
LimeSodaAsh_NCAR = csvread('LimeSodaAsh_NCAR.csv');
LimeSodaAsh_RES = csvread('LimeSodaAsh_RES.csv');
LimeSodaAsh_ETX = csvread('LimeSodaAsh_ETX.csv');

%Loading sulfate model outputs
Sulfate_AP = csvread('Sulfate_AP.csv');
Sulfate_EP = csvread('Sulfate_EP.csv');
Sulfate_GWP = csvread('Sulfate_GWP.csv');
Sulfate_ODP = csvread('Sulfate_ODP.csv');
Sulfate_POCP = csvread('Sulfate_POCP.csv');
Sulfate_PEU = csvread('Sulfate_PEU.csv');
Sulfate_CAR = csvread('Sulfate_CAR.csv');
Sulfate_NCAR = csvread('Sulfate_NCAR.csv');
Sulfate_RES = csvread('Sulfate_RES.csv');
Sulfate_ETX = csvread('Sulfate_ETX.csv');

%loading caustic soda model outputs
CausticSoda_AP = csvread('CausticSoda_AP.csv');
CausticSoda_EP = csvread('CausticSoda_EP.csv');
CausticSoda_GWP = csvread('CausticSoda_GWP.csv');
CausticSoda_ODP = csvread('CausticSoda_ODP.csv');
CausticSoda_POCP = csvread('CausticSoda_POCP.csv');
CausticSoda_PEU = csvread('CausticSoda_PEU.csv');
CausticSoda_CAR = csvread('CausticSoda_CAR.csv');
CausticSoda_NCAR = csvread('CausticSoda_NCAR.csv');
CausticSoda_RES = csvread('CausticSoda_RES.csv');
CausticSoda_ETX = csvread('CausticSoda_ETX.csv');

%Calculating the averages of the IEPSS model output
xlswrite('Output_Averages', mean(IEPSS_AP), 'Sheet 1', 'A3')
xlswrite('Output_Averages', mean(IEPSS_EP), 'Sheet 1', 'B3')
xlswrite('Output_Averages', mean(IEPSS_GWP), 'Sheet 1', 'C3')
xlswrite('Output_Averages', mean(IEPSS_ODP), 'Sheet 1', 'D3')
xlswrite('Output_Averages', mean(IEPSS_POCP), 'Sheet 1', 'E3')
xlswrite('Output_Averages', mean(IEPSS_PEU), 'Sheet 1', 'F3')
xlswrite('Output_Averages', mean(IEPSS_CAR), 'Sheet 1', 'G3')
xlswrite('Output_Averages', mean(IEPSS_NCAR), 'Sheet 1', 'H3')
xlswrite('Output_Averages', mean(IEPSS_RES), 'Sheet 1', 'I3')
xlswrite('Output_Averages', mean(IEPSS_ETX), 'Sheet 1', 'J3')
%mean(IEPSS_Mg_ends)

%Calculating the averages of the lime soda ash model output
xlswrite('Output_Averages', mean(LimeSodaAsh_AP), 'Sheet 1', 'A6')
xlswrite('Output_Averages', mean(LimeSodaAsh_EP), 'Sheet 1', 'B6')
xlswrite('Output_Averages', mean(LimeSodaAsh_GWP), 'Sheet 1', 'C6')
xlswrite('Output_Averages', mean(LimeSodaAsh_ODP), 'Sheet 1', 'D6')
xlswrite('Output_Averages', mean(LimeSodaAsh_POCP), 'Sheet 1', 'E6')
xlswrite('Output_Averages', mean(LimeSodaAsh_PEU), 'Sheet 1', 'F6')
xlswrite('Output_Averages', mean(LimeSodaAsh_CAR), 'Sheet 1', 'G6')
xlswrite('Output_Averages', mean(LimeSodaAsh_NCAR), 'Sheet 1', 'H6')
xlswrite('Output_Averages', mean(LimeSodaAsh_RES), 'Sheet 1', 'I6')
xlswrite('Output_Averages', mean(LimeSodaAsh_ETX), 'Sheet 1', 'J6')

%Calculating the averages of the sulfate model output
xlswrite('Output_Averages', mean(Sulfate_AP), 'Sheet 1', 'A9')
xlswrite('Output_Averages', mean(Sulfate_EP), 'Sheet 1', 'B9')
xlswrite('Output_Averages', mean(Sulfate_GWP), 'Sheet 1', 'C9')
xlswrite('Output_Averages', mean(Sulfate_ODP), 'Sheet 1', 'D9')
xlswrite('Output_Averages', mean(Sulfate_POCP), 'Sheet 1', 'E9')
xlswrite('Output_Averages', mean(Sulfate_PEU), 'Sheet 1', 'F9')
xlswrite('Output_Averages', mean(Sulfate_CAR), 'Sheet 1', 'G9')
xlswrite('Output_Averages', mean(Sulfate_NCAR), 'Sheet 1', 'H9')
xlswrite('Output_Averages', mean(Sulfate_RES), 'Sheet 1', 'I9')
xlswrite('Output_Averages', mean(Sulfate_ETX), 'Sheet 1', 'J9')

%Calculating the averages of the caustic soda model output
xlswrite('Output_Averages', mean(CausticSoda_AP), 'Sheet 1', 'A12')
xlswrite('Output_Averages', mean(CausticSoda_EP), 'Sheet 1', 'B12')
xlswrite('Output_Averages', mean(CausticSoda_GWP), 'Sheet 1', 'C12')
xlswrite('Output_Averages', mean(CausticSoda_ODP), 'Sheet 1', 'D12')
xlswrite('Output_Averages', mean(CausticSoda_POCP), 'Sheet 1', 'E12')
xlswrite('Output_Averages', mean(CausticSoda_PEU), 'Sheet 1', 'F12')
xlswrite('Output_Averages', mean(CausticSoda_CAR), 'Sheet 1', 'G12')
xlswrite('Output_Averages', mean(CausticSoda_NCAR), 'Sheet 1', 'H12')
xlswrite('Output_Averages', mean(CausticSoda_RES), 'Sheet 1', 'I12')
xlswrite('Output_Averages', mean(CausticSoda_ETX), 'Sheet 1', 'J12')
