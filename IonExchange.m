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

%% WATER CONSTITUENT FILES

%Loading the contaminant data (mol/L)
load('Ba_input')
load('Sr_input')
load('Ca_input')
load('Mg_input')
load('Na_input')

%=Milliequilavent conversion factors
Ba_meq=2/137.33;
Sr_meq=2/87.62;
Ca_meq=2/40.078;
Mg_meq=1/24.305;
Na_meq=1/22.990;

%Converting the contaminant data to milliequivalents
Ba_Inf = Ba_input * Ba_meq*100; 
Sr_Inf = Sr_input * Sr_meq*100;
Ca_Inf = Ca_input * Ca_meq*100;
Mg_Inf = Mg_input * Mg_meq*100;
Na_Inf = Na_input * Na_meq*100;

%Number of resin sites (meq/L) in each column
TR = 1000; 

%Equilibrium constant (K) value between barium and sodium (KBS), strontium
%and sodium (KSN), calcium and sodium (KCN), and magnesium and sodium (KMN)
%for PSS
KBN = 0.45; 
KSN = 0.32; 
KCN = 0.30;
KMN = 0.23; 

%Equilibrium constant (K) value between barium and sodium (KBS), strontium
%and sodium (KSN), calcium and sodium (KCN), and magnesium and sodium (KMN)
%for PAA
%KBN = 0.; 
%KSN = 0.; 
%KCN = 0.;
%KMN = 0.; 

%Contaminant limits (meq/L) (CHECK THESE)
Ca_limit = 5 * 10^(-4) * 2;
Mg_limit = 4.11 * 10^(-4) * 2;
Ba_limit = 1.46 * 10^(-5) * 2;
Sr_limit = 1.71 * 10^(-5) * 2;

%m = number of segments the column is divided into + 1 (for initial
%conditions). n = initial estimate for the number of bed volumes treated
%(code stops on breakthrough of contaminant at a given level set above). 
%l = number of data points for each contaminant.
m = 10; 
n = 100*m;
l = 3;
bed_volumes = zeros(l,1);

%% SOLVER FUNCTION
for k=1:l
    
%Preallocating resin and water variables and/or filling them with zeroes 
    RESIN = zeros(n,m);
    WATER = zeros(n,m);
    
    for i=1:n
        for j=1:m
            
%Initial condition of virgin resin (preloaded with Na)                                   
            RESIN.Na(1,1:m) = TR;
            RESIN.Ba(1,1:m) = 0;
            RESIN.Sr(1,1:m) = 0;
            RESIN.Ca(1,1:m) = 0;
            RESIN.Mg(1,1:m) = 0;
            
%Initial condition of water 
            WATER.Na(1:n,1)=Na_Inf(k);
            WATER.Ba(1:n,1)=Ba_Inf(k);
            WATER.Sr(1:n,1)=Sr_Inf(k);
            WATER.Ca(1:n,1)=Ca_Inf(k);
            WATER.Mg(1:n,1)=Mg_Inf(k);
                       
%Storing the matrix data into more convenient (shorter) variables
            RNa = RESIN.Na(i,j);
            RBa = RESIN.Ba(i,j);
            RSr = RESIN.Sr(i,j);
            RCa = RESIN.Ca(i,j);
            RMg = RESIN.Mg(i,j);        
            WNa = WATER.Na(i,j);
            WBa = WATER.Ba(i,j);
            WSr = WATER.Sr(i,j);
            WCa = WATER.Ca(i,j);
            WMg = WATER.Mg(i,j);
       
            set0=[RNa;RBa;RSr;RCa;RMg;WNa;WBa;WSr;WCa;WMg];

            TNa = RNa + WNa;
            TBa = RBa + WBa;
            TSr = RSr + WSr;
            TCa = RCa + WCa; 
            TMg = RMg + WMg;

            options = optimset('Display','off','TolFun', 1.0e-4, 'TolX',1.0e-4);
            f = @(dummy)msolve(dummy,TR,KBN,KSN,KCN,KMN,TNa,TBa,TSr,TCa,TMg);

            [set] = fsolve(f,set0,options);

            RESIN.Na(i+1,j) = set(1);
            RESIN.Ba(i+1,j) = set(2);
            RESIN.Sr(i+1,j) = set(3);
            RESIN.Ca(i+1,j) = set(4);
            RESIN.Mg(i+1,j) = set(5);
            WATER.Na(i,j+1) = set(6);
            WATER.Ba(i,j+1) = set(7);
            WATER.Sr(i,j+1) = set(8);
            WATER.Ca(i,j+1) = set(9);
            WATER.Mg(i,j+1) = set(10);

        end
        
%Ends the modeling if one of the contaminants is above the set limit
             if WBa > Ba_limit || WSr > Sr_limit || WCa > Ca_limit || WMg > Mg_limit
                  break
             end 
            
%Stores the number of bed volumes completed             
            bed_volumes(k,1) = floor(i/m)
    end

end
                                                                 
%Saves the bed volumes, rounded down to the nearest whole number, to a file
xlswrite('IonExchange_BedVolumes.xlsx', bed_volumes);