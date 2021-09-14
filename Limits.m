%Calculates the Ca, Ba, and Sr limits based on RO treatment goal and SO4

load('SO4_input'); %as mol/l

%Converting SO4_input to mg/l
SO4_input = SO4_input * 96.056 * 1000;

%Constants
mass = SO4_input; %Mass
vol = 1; %Volume
sel = 0.95; %SO4 selectivity
rec = 0.5; %Recovery
Ksp_Ca = 6.91831*10^(-5); %Ca selectivity
Ksp_Ba = 1.0842*10^(-10); %Ba selectivity
Ksp_Sr = 3.44*10^(-7); %Sr selectivity 
l = 30000; %number of data points
Ca_ends = zeros(l,1);
Ba_ends = zeros(l,1);
Sr_ends = zeros(l,1);
SO4 = zeros(l,1);

%Main loop
for i=1:l
    
    perm = (1-sel)*mass(i);

    rej = (mass(i)-perm)/((1-rec)*vol); %mg/l

    rej = rej / (96*1000); %mol/l

    Ca_sat = Ksp_Ca/rej; %mol/l
    Ca_ends(i) = Ca_sat*(1-rec)/sel; %mol/l
    Ca_ends(i) = 2*1000*Ca_ends(i); %meq/l
    
    Ba_sat = Ksp_Ba/rej; %mol/l
    Ba_ends(i) = Ba_sat*(1-rec)/sel; %mol/l
    Ba_ends(i) = 2*1000*Ba_ends(i); %meq/l
        
    Sr_sat = Ksp_Sr/rej; %mol/l
    Sr_ends(i) = Sr_sat*(1-rec)/sel; %mol/l
    Sr_ends(i) = 2*1000*Sr_ends(i); %meq/l

end

%Saving limits to output file
csvwrite('Ca_ends', Ca_ends);
csvwrite('Ba_ends', Ba_ends);
csvwrite('Sr_ends', Sr_ends);