%This script finds all influent waters that don't need treatment

load('Ca_input')
load('Mg_input')
load('Ba_input')
load('Sr_input')
load('Alk_input')
load('Mg_ends')
load('Ca_ends')
load('Sr_ends')
load('Ba_ends')

%n is the number of data points in the data set, count is the number of
%waters not needing treatment, m is the number of data points that will be
%in the new influent water data set (which all require some treatment), l
%keeps track of the new data set index. Mg is not used for this because
%only ~2% of the time Mg is the only required contaminant and the limit
%calculation of Mg for 26k data points would take ~15 hrs vie IEPSS_blended
%model
n = 30000;
count = 0;
m = 10000;
l = 1;

%Creating non-zero influent water vectors and their corresponding limits
Ca_input_nz = zeros(m,1);
Mg_input_nz = zeros(m,1);
Ba_input_nz = zeros(m,1);
Sr_input_nz = zeros(m,1);
Alk_input_nz = zeros(m,1);
Ca_ends_nz = zeros(m,1);
Mg_ends_nz = zeros(m,1);
Ba_ends_nz = zeros(m,1);
Sr_ends_nz = zeros(m,1);

%This section runs through the data and displays the position and running
%total of waters not needing treatment
for i = 1:n
    if Ca_input(i) < Ca_ends(i) && Ba_input(i) < Ba_ends(i) && Sr_input(i) < Sr_ends(i)
        %disp(i)
        %count = count + 1;
        %disp(count)
    end
end


%This section finds only waters that require some treatment and adds them
%into a new variable along with their corresponding limits
for i = 1:n
    if Ca_input(i) > Ca_ends(i) || Ba_input(i) > Ba_ends(i) || Sr_input(i) > Sr_ends(i)
        Ca_input_nz(l) = Ca_input(i);
        Mg_input_nz(l) = Mg_input(i);
        Ba_input_nz(l) = Ba_input(i);
        Sr_input_nz(l) = Sr_input(i);
        Alk_input_nz(l) = Alk_input(i);
        Ca_ends_nz(l) = Ca_ends(i);
        Ba_ends_nz(l) = Ba_ends(i);
        Sr_ends_nz(l) = Sr_ends(i);
        l = l+1    
    end
    
    if l > m
        break;
    end
   
end

%This script still needs to keep track of the water sets that don't need
%treatment and remove them from the original and create a new set of only
%waters needing treatment