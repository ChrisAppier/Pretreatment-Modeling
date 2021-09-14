%% ION EXCHANGE MODEL - ONLY ACCOUNTS FOR OPERATION, NO RENGERATION OF THE COLUMN

%Code solves for the equilibrium of cations in water simultaneous 

%Code divides an ion exchange column into a set number of segments (m).
%These segments are defined by their capacity as meq/L. A meq is charged
%based measurement of the cations in water.

%The incoming water is divided into similar volume size as the segments in
%the column and are known as packets of water (n). Each packet of water is
%assumed to be in instantaneous equilibrium with the column. The
%equilibrium of the cations with the column is determined by the equilbrium
%coefficients (K). Five cations are studied in this experiment - sodium,
%barium, strontium, calcium, and magnesium.

tic
%warning('off', 'all')

%% WATER CONSTITUENT FILES

%Loading the contaminant data (mol/L)
load('Ba_input')
load('Sr_input')
load('Ca_input')
load('Mg_input')
load('Na_input')
load('Ca_ends')
load('Sr_ends')
load('Ba_ends')

%=Molarity to Milliequilavent/L conversion factors
Ba_meq = 2 * 1000;
Sr_meq = 2 * 1000;
Ca_meq = 2 * 1000;
Mg_meq = 2 * 1000;
Na_meq = 1 * 1000;

%Converting the influent data to milliequivalents
Ba_Inf = Ba_input * Ba_meq;
Sr_Inf = Sr_input * Sr_meq;
Ca_Inf = Ca_input * Ca_meq;
Mg_Inf = Mg_input * Mg_meq;
Na_Inf = Na_input * Na_meq;

%Number of resin sites (meq/L) in each column
TR = 4325; 

%Equilibrium constant (K) value between barium and sodium (KBS), strontium
%and sodium (KSN), calcium and sodium (KCN), and magnesium and sodium (KMN)
%for PSS
KBN = 5.41; %values obtained from DuPont 
KSN = 2.62; 
KCN = 1.89;
KMN = 1.48;

%Contaminant limits (meq/L)
Ca_limit = Ca_ends;
Mg_limit = 10000000; %no limit
Ba_limit = Ba_ends;
Sr_limit = Sr_ends;

%Defining the LCA data values
PSS_AP = 45.38*0.00706 + 0.513*0.0117; %Acidification potential (kg SO2 eq per kg PSS)
PSS_EP = 45.38*0.000729 + 0.513*0.0012; %Eutrophication potential (kg N eq per kg PSS)
PSS_GWP = 45.38*0.098 + 0.513*3.59; %Global warming potential (kg CO2 eq per kg PSS)
PSS_ODP = 45.38*4.59*10^(-8) + 0.513*2.05*10^(-8); %Ozone depletion potential (kg CFC-11 eq per kg PSS)
PSS_POCP = 45.38*0.0187 + 0.513*0.146; %Photochemical ozone creation potential (kg O3 eq per kg PSS)
PSS_PEU = 45.38*0.979 + 0.513*12.1; %Primary energy use (MJ surplus per kg PSS)
PSS_CAR = 45.38*1.48*10^(-8) + 0.513*9.15*10^(-8); %Carcinogenics (CTUh per kg PSS)
PSS_NCAR = 45.38*0.000000124 + 0.513*0.000000107; %Non-carcinogenics (CTUh per kg PSS)
PSS_RES = 45.38*0.000496 + 0.513*0.0000904; %Respiratory effects (kg PM2.5 eq per kg PSS)
PSS_ETX = 45.38*3.03 + 0.513*7.84; %Ecotoxicity (CTUe per kg PSS)

%m = number of segments the column is divided into + 1 (for initial
%conditions). n = initial estimate for the number of bed volumes treated
%(code stops on breakthrough of contaminant at a given level set above). 
%l = number of data points for each contaminant.
m = 20;
n = m*50;
l = 10000;
bed_volumes = zeros(l,1);
Mg_ends = zeros(l,1);
AP = zeros(l,1);
EP = zeros(l,1);
GWP = zeros(l,1);
ODP = zeros(l,1);
POCP = zeros(l,1);
PEU = zeros(l,1);
CAR = zeros(l,1);
NCAR = zeros(l,1);
RES = zeros(l,1);
ETX = zeros(l,1);

resin_used = zeros(l,1);
res_den = 801; %g/l Sigma Aldrich

%% SOLVER FUNCTION
parfor k=1:l
    disp(k)
%Preallocating/resetting resin and water variables and filling them with zeroes 
    WATERNa = zeros(n,m);
    WATERBa = zeros(n,m);
    WATERSr = zeros(n,m);
    WATERCa = zeros(n,m);
    WATERMg = zeros(n,m);
    RESINNa = zeros(n,m);
    RESINBa = zeros(n,m);
    RESINSr = zeros(n,m);
    RESINCa = zeros(n,m);
    RESINMg = zeros(n,m);
    Ba_avg = zeros(1,n);
    Sr_avg = zeros(1,n);
    Ca_avg = zeros(1,n);
    Mg_avg = zeros(1,n);
            
    for i=1:n

        if i == n
            disp('bruh')
        end
        
        for j=1:m
%Initial condition of virgin resin (preloaded with Na)                                               
            RESINNa(1,1:m) = TR;
            RESINBa(1,1:m) = 0;
            RESINSr(1,1:m) = 0;
            RESINCa(1,1:m) = 0;
            RESINMg(1,1:m) = 0;
            
%Initial condition of water 
            WATERNa(1:n,1)=Na_Inf(k);
            WATERBa(1:n,1)=Ba_Inf(k);
            WATERSr(1:n,1)=Sr_Inf(k);
            WATERCa(1:n,1)=Ca_Inf(k);
            WATERMg(1:n,1)=Mg_Inf(k);
                       
%Storing the matrix data into more convenient (shorter) variables
            RNa = RESINNa(i,j);
            RBa = RESINBa(i,j);
            RSr = RESINSr(i,j);
            RCa = RESINCa(i,j);
            RMg = RESINMg(i,j);        
            WNa = WATERNa(i,j);
            WBa = WATERBa(i,j);
            WSr = WATERSr(i,j);
            WCa = WATERCa(i,j);
            WMg = WATERMg(i,j);
       
%Solving the system of linear equations
            set0=[RNa;RBa;RSr;RCa;RMg;WNa;WBa;WSr;WCa;WMg];

            TNa = RNa + WNa;
            TBa = RBa + WBa;
            TSr = RSr + WSr;
            TCa = RCa + WCa; 
            TMg = RMg + WMg;

            options = optimset('Display','off','TolFun', 1.0e-4, 'TolX',1.0e-4);
            
            f = @(dummy)msolve(dummy,TR,KBN,KSN,KCN,KMN,TNa,TBa,TSr,TCa,TMg);

            [set] = fsolve(f,set0,options);

%Assigning the new values calculated above            
            RESINNa(i+1,j) = set(1);
            RESINBa(i+1,j) = set(2);
            RESINCa(i+1,j) = set(4);
            RESINMg(i+1,j) = set(5);
            WATERNa(i,j+1) = set(6);
            WATERBa(i,j+1) = set(7);
            WATERSr(i,j+1) = set(8);
            WATERCa(i,j+1) = set(9);
            WATERMg(i,j+1) = set(10);
            
        end
        
%Calculating the average contaminant levels of all of the treated water        
            Ba_avg_temp = mean(WATERBa);
            Ba_avg(i) = Ba_avg_temp(1,m) * (n/i); %(n/i) accounts for pre-filled zeros
            Ca_avg_temp = mean(WATERCa);
            Ca_avg(i) = Ca_avg_temp(1,m) * (n/i);
            Sr_avg_temp = mean(WATERSr);
            Sr_avg(i) = Sr_avg_temp(1,m) * (n/i);
            Mg_avg_temp = mean(WATERMg);
            Mg_avg(i) = Mg_avg_temp(1,m) * (n/i);
            

%Ends the modeling if one of the contaminants is above the set limit and
%announces the concentration
            if Ba_avg(i) > Ba_limit(k,1) || Sr_avg(i) > Sr_limit(k,1) || Ca_avg(i) > Ca_limit(k,1)
                 break
            end 
            
%Stores the number of bed volumes completed and Mg end          
            bed_volumes(k,1) = i/m;
            Mg_ends(k,1) = Mg_avg(i);
    end

%Calculates amount of resin used per volume water (kg/m^3)
    resin_used(k,1) = res_den * 1000 / bed_volumes(k,1);

%Calculates emissions impacts
    AP(k,1) = resin_used(k,1)*PSS_AP; %Acidification potential
    EP(k,1) = resin_used(k,1)*PSS_EP; %Eutrophication potential
    GWP(k,1) = resin_used(k,1)*PSS_GWP; %Global warming potential
    ODP(k,1) = resin_used(k,1)*PSS_ODP; %Ozone depletion potential
    POCP(k,1) = resin_used(k,1)*PSS_POCP; %Photochemical ozone creation potential
    PEU(k,1) = resin_used(k,1)*PSS_PEU; %Primary energy use
    CAR(k,1) = resin_used(k,1)*PSS_CAR; %Carcinogenics
    NCAR(k,1) = resin_used(k,1)*PSS_CAR; %Non-carcinogenics
    RES(k,1) = resin_used(k,1)*PSS_CAR; %Respiratory effects
    ETX(k,1) = resin_used(k,1)*PSS_CAR; %Ecotoxicity
end
%% Debug

%Outputting final answer to command window
%BV = bed_volumes;

%Plotting initial water packet pollutant concentrations moving through the
%column number set below
%col_num = 1; %Set equal to mBV to see the last column
%for i = 1:m
    %Ba(i,1) = WATERBa(col_num,i);
    %Sr(i,1) = WATERSr(col_num,i);
    %Ca(i,1) = WATERCa(col_num,i);
    %Mg(i,1) = WATERMg(col_num,i);
    %Na(i,1) = WATERNa(col_num,i);
%end

%figure
%subplot(1,5,1); plot(Ba, 'Linewidth', 2.0); title('Ba (aq)'); xlabel('Segments'); ylabel('Concentration [meq/L]');
%subplot(1,5,2); plot(Sr, 'Linewidth', 2.0); title('Sr (aq)'); xlabel('Segments'); ylabel('Concentration [meq/L]');
%subplot(1,5,3); plot(Ca, 'Linewidth', 2.0); title('Ca (aq)'); xlabel('Segments'); ylabel('Concentration [meq/L]');
%subplot(1,5,4); plot(Mg, 'Linewidth', 2.0); title('Mg (aq)'); xlabel('Segments'); ylabel('Concentration [meq/L]');
%subplot(1,5,5); plot(Na, 'Linewidth', 2.0); title('Na (aq)'); xlabel('Segments'); ylabel('Concentration [meq/L]');

%Plotting the water concentrations after each column and the concentration
%limits
%for i = 1:n
    %Ba_final(i) = WATERBa(i,m);
    %Sr_final(i) = WATERSr(i,m);
    %Ca_final(i) = WATERCa(i,m);
    %Mg_final(i) = WATERMg(i,m);
    %Na_final(i) = WATERNa(i,m);
    
    %Ba_lim(i,1) = Ba_limit;
    %Ca_lim(i,1) = Ca_limit;
    %Sr_lim(i,1) = Sr_limit;
    %Mg_lim(i,1) = Mg_limit;
    
    %column(i) = i;
%end

%figure
%subplot(1,5,1); plot(column, Ba_final, column, Ba_lim, '--', 'Linewidth', 2.0); title('Ba (aq)'); xlim([0 m*BV+1]);  xlabel('BV'); ylabel('Concentration [meq/L]');
%subplot(1,5,2); plot(column, Sr_final, column, Sr_lim, '--', 'Linewidth', 2.0); title('Sr (aq)'); xlim([0 m*BV+1]);xlabel('BV'); ylabel('Concentration [meq/L]');
%subplot(1,5,3); plot(column, Ca_final, column, Ca_lim, '--', 'Linewidth', 2.0); title('Ca (aq)'); xlim([0 m*BV+1]);xlabel('BV'); ylabel('Concentration [meq/L]');
%subplot(1,5,4); plot(column, Mg_final, column, Mg_lim, '--', 'Linewidth', 2.0); title('Mg (aq)'); xlim([0 m*BV+1]);xlabel('BV'); ylabel('Concentration [meq/L]');
%subplot(1,5,5); plot(Na_final, 'Linewidth', 2.0); title('Na (aq)', 'Linewidth', 2.0); xlim([0 mBV+1]);xlabel('BV * m'); ylabel('Concentration [meq/L]');

%Plotting the average pollutant concentration in all treated water
%figure
%subplot(1,4,1); plot(column,Ba_avg, column, Ba_lim, '--', 'Linewidth', 2.0); title('Avg Ba (aq)'); xlim([0 m*BV]);
%subplot(1,4,2); plot(column,Sr_avg, column, Sr_lim, '--', 'Linewidth', 2.0); title('Avg Sr (aq)'); xlim([0 m*BV]);
%subplot(1,4,3); plot(column,Ca_avg, column, Ca_lim, '--', 'Linewidth', 2.0); title('Avg Ca (aq)'); xlim([0 m*BV]);
%subplot(1,4,4); plot(column,Mg_avg, column, Mg_lim, '--', 'Linewidth', 2.0); title('Avg Mg (aq)'); xlim([0 m*BV]);

%Saving water and resin matrices to file for debug
%csvwrite('Ca_Water_Matrix.csv', WATERCa);
%csvwrite('Ba_Water_Matrix.csv', WATERBa);
%csvwrite('Na_Water_Matrix.csv', WATERNa);
%csvwrite('Sr_Water_Matrix.csv', WATERSr);
%csvwrite('Mg_Water_Matrix.csv', WATERMg);
%csvwrite('Ca_Resin_Matrix.csv', RESINCa);
%csvwrite('Ba_Resin_Matrix.csv', RESINBa);
%csvwrite('Na_Resin_Matrix.csv', RESINNa);
%csvwrite('Sr_Resin_Matrix.csv', RESINSr);
%csvwrite('Mg_Resin_Matrix.csv', RESINMg);


%% Output 
%Saves the bed volumes, rounded down to the nearest whole number, to a file
%Also saves Mg_limit for use in other models
csvwrite('IEPSS_BedVolumes.csv', bed_volumes);
csvwrite('Mg_ends', Mg_ends);
csvwrite('IEPSS_AP.csv', AP);
csvwrite('IEPSS_EP.csv', EP);
csvwrite('IEPSS_GWP.csv', GWP);
csvwrite('IEPSS_ODP.csv', ODP);
csvwrite('IEPSS_POCP.csv', POCP);
csvwrite('IEPSS_PEU.csv', PEU);
csvwrite('IEPSS_CAR.csv', CAR);
csvwrite('IEPSS_NCAR.csv', NCAR);
csvwrite('IEPSS_RES.csv', RES);
csvwrite('IEPSS_ETX.csv', ETX);

%Writes average (median) environmental impact factors to screen
disp(median(AP))
disp(median(EP))
disp(median(GWP))
disp(median(ODP))
disp(median(POCP))
disp(median(PEU))
disp(median(CAR))
disp(median(NCAR))
disp(median(RES))
disp(median(ETX))

toc