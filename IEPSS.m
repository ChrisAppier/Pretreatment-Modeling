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

clc
clear
tic
%warning('off', 'all')

%% WATER CONSTITUENT FILES

%Loading the contaminant data (mol/L)
load('Ba_input')
load('Sr_input')
load('Ca_input')
load('Mg_input')
load('Na_input')

%=Molarity to Milliequilavent/L conversion factors
Ba_meq = 2 * 1000;
Sr_meq = 2 * 1000;
Ca_meq = 2 * 1000;
Mg_meq = 1 * 1000;
Na_meq = 1 * 1000;

%Converting the influent data to milliequivalents
Ba_Inf = Ba_input * Ba_meq;
Sr_Inf = Sr_input * Sr_meq;
Ca_Inf = Ca_input * Ca_meq;
Mg_Inf = Mg_input * Mg_meq;
Na_Inf = Na_input * Na_meq;

%Number of resin sites (meq/L) in each column
TR = 385; 

%Equilibrium constant (K) value between barium and sodium (KBS), strontium
%and sodium (KSN), calcium and sodium (KCN), and magnesium and sodium (KMN)
%for PSS
KBN = 0.45; 
KSN = 0.32; 
KCN = 0.30;
KMN = 0.23; 

%Contaminant limits (meq/L) (Condition 1)
Ca_limit = 8.234 * 10^(-3) * Ca_meq;
Mg_limit = 6.418 * 10^(-3) * Mg_meq;
Ba_limit = 1.566 * 10^(-4) * Ba_meq;
Sr_limit = 9.130 * 10^(-5) * Sr_meq;

%Contaminant limits (meq/L) (Condition 2)
%Ca_limit = 1.153 * 10^(-1) * Ca_meq;
%Mg_limit = 5.925 * 10^(-2) * Mg_meq;
%Ba_limit = 1.000 * 10^(14) * Ba_meq;
%Sr_limit = 1.664 * 10^(-2) * Sr_meq;

%Contaminant limits (meq/L) (Condition 3)
%Ca_limit = 1.347 * 10^(-2) * Ca_meq;
%Mg_limit = 8.825 * 10^(-3) * Mg_meq;
%Ba_limit = 1.000 * 10^(14) * Ba_meq;
%Sr_limit = 1.826 * 10^(-5) * Sr_meq;

%Debug
Ca_change = Ca_Inf - Ca_limit;
Ba_change = Ba_Inf - Ba_limit;
Mg_change = Mg_Inf - Mg_limit;
Sr_change = Sr_Inf - Sr_limit;

%m = number of segments the column is divided into + 1 (for initial
%conditions). n = initial estimate for the number of bed volumes treated
%(code stops on breakthrough of contaminant at a given level set above). 
%l = number of data points for each contaminant.
m = 100;
n = m*2;
l = 1;
bed_volumes = zeros(l,1);

%% SOLVER FUNCTION
for k=l:l
    
%Preallocating resin and water variables and/or filling them with zeroes 
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
    
    for i=1:n
        i = i
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
       
            set0=[RNa;RBa;RSr;RCa;RMg;WNa;WBa;WSr;WCa;WMg];

            TNa = RNa + WNa;
            TBa = RBa + WBa;
            TSr = RSr + WSr;
            TCa = RCa + WCa; 
            TMg = RMg + WMg;

            options = optimset('Display','off','TolFun', 1.0e-4, 'TolX',1.0e-4);
            f = @(dummy)msolve(dummy,TR,KBN,KSN,KCN,KMN,TNa,TBa,TSr,TCa,TMg);

            [set] = fsolve(f,set0,options);

            RESINNa(i+1,j) = set(1);
            RESINBa(i+1,j) = set(2);
            RESINSr(i+1,j) = set(3);
            RESINCa(i+1,j) = set(4);
            RESINMg(i+1,j) = set(5);
            WATERNa(i,j+1) = set(6);
            WATERBa(i,j+1) = set(7);
            WATERSr(i,j+1) = set(8);
            WATERCa(i,j+1) = set(9);
            WATERMg(i,j+1) = set(10);
            
        end
        
%Ends the modeling if one of the contaminants is above the set limit
             %if WBa > Ba_limit || WSr > Sr_limit || WCa > Ca_limit || WMg > Mg_limit
                  %break
             %end 
             
             if WBa > Ba_limit
                 disp('Ba =')
                 disp(WBa)
                 break
             end
             
             if WCa > Ca_limit
                 disp('Ca =')
                 disp(WCa)
                 break
             end
            
             if WSr > Sr_limit
                 disp('Sr =')
                 disp(WSr)
                 break
             end
             
             if WMg > Mg_limit
                 disp('Ba =')
                 disp(WMg)
                 break
             end
             
%Stores the number of bed volumes completed             
            %bed_volumes(k,1) = floor(i/m);
            bed_volumes(k,1) = i;         
    end

end

%% Output

%Outputting final answer to command window
mBV = bed_volumes

%Plotting initial water packet pollutant concentrations moving through the first column
for i = 1:m
    Ba(i,1) = WATERBa(l,i);
    Sr(i,1) = WATERSr(l,i);
    Ca(i,1) = WATERCa(l,i);
    Mg(i,1) = WATERMg(l,i);
    Na(i,1) = WATERNa(l,i);
end

figure
subplot(1,5,1); plot(Ba, 'Linewidth', 2.0); title('Ba (aq)'); xlabel('Segments'); ylabel('Concentration [meq/L]');
subplot(1,5,2); plot(Sr, 'Linewidth', 2.0); title('Sr (aq)'); xlabel('Segments'); ylabel('Concentration [meq/L]');
subplot(1,5,3); plot(Ca, 'Linewidth', 2.0); title('Ca (aq)'); xlabel('Segments'); ylabel('Concentration [meq/L]');
subplot(1,5,4); plot(Mg, 'Linewidth', 2.0); title('Mg (aq)'); xlabel('Segments'); ylabel('Concentration [meq/L]');
subplot(1,5,5); plot(Na, 'Linewidth', 2.0); title('Na (aq)'); xlabel('Segments'); ylabel('Concentration [meq/L]');

%Plotting the final water concentrations after each column
for i = 1:n
    Ba_final(i) = WATERBa(i,m);
    Sr_final(i) = WATERSr(i,m);
    Ca_final(i) = WATERCa(i,m);
    Mg_final(i) = WATERMg(i,m);
    Na_final(i) = WATERNa(i,m);
    
    Ba_lim(i,1) = Ba_limit;
    Ca_lim(i,1) = Ca_limit;
    Sr_lim(i,1) = Sr_limit;
    Mg_lim(i,1) = Mg_limit;
    
    column(i) = i;
end

figure
subplot(1,5,1); plot(column, Ba_final, column, Ba_lim, '--', 'Linewidth', 2.0); title('Ba (aq)'); ylim([0 1.1*max(Ba_final)]); xlim([0 mBV(l,1)+1]);  xlabel('BV * m'); ylabel('Concentration [meq/L]');
subplot(1,5,2); plot(column, Sr_final, column, Sr_lim, '--', 'Linewidth', 2.0); title('Sr (aq)'); ylim([0 1.1*max(Sr_final)]); xlim([0 mBV(l,1)+1]);xlabel('BV * m'); ylabel('Concentration [meq/L]');
subplot(1,5,3); plot(column, Ca_final, column, Ca_lim, '--', 'Linewidth', 2.0); title('Ca (aq)'); ylim([0 1.1*max(Ca_final)]); xlim([0 mBV(l,1)+1]);xlabel('BV * m'); ylabel('Concentration [meq/L]');
subplot(1,5,4); plot(column, Mg_final, column, Mg_lim, '--', 'Linewidth', 2.0); title('Mg (aq)'); ylim([0 1.1*max(Mg_final)]); xlim([0 mBV(l,1)+1]);xlabel('BV * m'); ylabel('Concentration [meq/L]');
subplot(1,5,5); plot(Na_final, 'Linewidth', 2.0); title('Na (aq)', 'Linewidth', 2.0); xlim([0 mBV(l,1)+1]);xlabel('BV * m'); ylabel('Concentration [meq/L]');

%Saving water matrices to file for debug
%csvwrite('Ca_Water_Matrix.csv', WATERCa);
                                                                 
%Saves the bed volumes, rounded down to the nearest whole number, to a file
csvwrite('IEPSS_BedVolumes.csv', bed_volumes);

toc