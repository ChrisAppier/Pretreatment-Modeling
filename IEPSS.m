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
warning('off', 'all')

%% WATER CONSTITUENT FILES

%Loading the contaminant data (mol/L)
load('Ba_input')
load('Sr_input')
load('Ca_input')
load('Mg_input')
load('Na_input')

%=Milliequilavent conversion factors
Ba_meq = 2;
Sr_meq = 2;
Ca_meq = 2;
Mg_meq = 1;
Na_meq = 1;

%Converting the contaminant data to milliequivalents
Ba_Inf = Ba_input * Ba_meq; 
Sr_Inf = Sr_input * Sr_meq;
Ca_Inf = Ca_input * Ca_meq;
Mg_Inf = Mg_input * Mg_meq;
Na_Inf = Na_input * Na_meq;

%Number of resin sites (meq/L) in each column
TR = 208; 

%Equilibrium constant (K) value between barium and sodium (KBS), strontium
%and sodium (KSN), calcium and sodium (KCN), and magnesium and sodium (KMN)
%for PSS
KBN = 0.45; 
KSN = 0.32; 
KCN = 0.30;
KMN = 0.23; 

%Contaminant limits (meq/L)
Ca_limit = 8.234 * 10^(-3) * 2;
Mg_limit = 6.418 * 10^(-3) * 2;
Ba_limit = 1.566 * 10^(-4) * 2;
Sr_limit = 9.130 * 10^(-5) * 2;

%m = number of segments the column is divided into + 1 (for initial
%conditions). n = initial estimate for the number of bed volumes treated
%(code stops on breakthrough of contaminant at a given level set above). 
%l = number of data points for each contaminant.
m = 100;
n = 150*m;
l = 1;
bed_volumes = zeros(l,1);

%% SOLVER FUNCTION
parfor k=1:l
    
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
             if WBa > Ba_limit || WSr > Sr_limit || WCa > Ca_limit || WMg > Mg_limit
                  break
             end 
            
%Stores the number of bed volumes completed             
            bed_volumes(k,1) = floor(i/m);
    end

end

%% Output

%Outputting final answer to command window
BV = bed_volumes
                                                                 
%Saves the bed volumes, rounded down to the nearest whole number, to a file
csvwrite('IEPSS_BedVolumes.csv', bed_volumes);

toc