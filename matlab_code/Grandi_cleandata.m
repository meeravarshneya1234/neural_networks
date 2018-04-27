%% Remove APs with EADs before the perturbation 
load('Grandi_1000rawpop.mat')

x = datatable.ctl.times;
y = datatable.ctl.V;

figure
hold on
cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
title('Original')

x = datatable.exp.times;
y = datatable.exp.V;

figure
hold on
cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
title('Original - after perturbation')

t_cutoff = 4;
flag = 0;
[APfails,EADs] = cleandata(datatable.ctl.APDs,datatable.ctl.times,datatable.ctl.V,t_cutoff,flag);
[number_of_failed,ind_failed] = find(APfails ==1); % number of failed to repolarize
[number_of_EADs,ind_EADs] = find(EADs==1); % number of EADs

clean_datatable = [];

clean_datatable.times(:,1) = datatable.ctl.times(~(EADs' + APfails'));
clean_datatable.V(:,1) = datatable.ctl.V(~(EADs' + APfails'));
clean_datatable.states(:,1) = datatable.ctl.states(~(EADs' + APfails'));
clean_datatable.APDs(:,1) = datatable.ctl.APDs(~(EADs' + APfails'));
clean_datatable.scaling = datatable.ctl.scalings(~(EADs' + APfails'),:);

clean_datatable.times(:,2) = datatable.exp.times(~(EADs' + APfails'));
clean_datatable.V(:,2) = datatable.exp.V(~(EADs' + APfails'));
clean_datatable.states(:,2) = datatable.exp.states(~(EADs' + APfails'));
clean_datatable.APDs(:,2) = datatable.exp.APDs(~(EADs' + APfails'));
clean_datatable.scaling = datatable.exp.scalings(~(EADs' + APfails'),:);

x = clean_datatable.times(:,1);
y = clean_datatable.V(:,1);

figure
hold on
cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
title('Cleaned Data - before perturbation')

x2 = clean_datatable.times(:,2);
y2 = clean_datatable.V(:,2);

figure
hold on
cellfun(@(x,y) plot(x,y,'linewidth',2),x2,y2)
title('Cleaned Data- after perturbation')

%% Separate APs based on whether EAD occured after the perturbation 
[APfails,EADs] = cleandata(clean_datatable.APDs(:,2),clean_datatable.times(:,2),clean_datatable.V(:,2),1);
%large_APDs = (clean_datatable.APDs(:,2)>700); % specifcally made for this case, AP had an EAD but was not detected 
% by cleandata function

with_EADs = [];
no_EADs = [];
for aF = 1:2
    I = ~(EADs' + APfails'); %no EAD 
    no_EADs.times(:,aF) = clean_datatable.times(I,aF);
    no_EADs.V(:,aF) = clean_datatable.V(I,aF);  
    no_EADs.states(:,aF) = clean_datatable.states(I,aF);  
    no_EADs.APDs(:,aF) = clean_datatable.APDs(I,aF); 
    
    with_EADs.times(:,aF) = clean_datatable.times(~I,aF);
    with_EADs.V(:,aF) = clean_datatable.V(~I,aF);  
    with_EADs.states(:,aF) = clean_datatable.states(~I,aF);  
    with_EADs.APDs(:,aF) = clean_datatable.APDs(~I,aF); 
        
end

% for aF = 1:2
%     I = ~(EADs' + APfails'+large_APDs);
%     no_EADs.times(:,aF) = clean_datatable.times(~(EADs' + APfails'+large_APDs),aF);
%     no_EADs.V(:,aF) = clean_datatable.V(~(EADs' + APfails'+large_APDs),aF);  
%     no_EADs.states(:,aF) = clean_datatable.states(~(EADs' + APfails'+large_APDs),aF);  
%     no_EADs.APDs(:,aF) = clean_datatable.APDs(~(EADs' + APfails'+large_APDs),aF); 
%     
%     with_EADs.times(:,aF) = clean_datatable.times(~I,aF);
%     with_EADs.V(:,aF) = clean_datatable.V(~I,aF);  
%     with_EADs.states(:,aF) = clean_datatable.states(~I,aF);  
%     with_EADs.APDs(:,aF) = clean_datatable.APDs(~I,aF); 
% end

x3 = no_EADs.times(:,2);
y3 = no_EADs.V(:,2);
x4 = with_EADs.times(:,2);
y4 = with_EADs.V(:,2);

figure
hold on
cellfun(@(x,y) plot(x,y,'linewidth',2),x3,y3)
title(['After perturbation APs with no EADs = ' num2str(length(x3))])

figure
hold on
cellfun(@(x,y) plot(x,y,'linewidth',2),x4,y4)
title(['After perturbation APs formed EADs = ' num2str(length(x4))])

%% 
NN_train_y = EADs + APfails; 
NN_train_y = NN_train_y';

t0 = 0:1:999;
ts = clean_datatable.times(:,1);
ys = clean_datatable.V(:,1);

for i = 1:length(ts) 
    y0(:,i) = interp1(ts{i},ys{i},t0);
end 
y0 = y0'; 
y0(:,end+1) = NN_train_y;

fid = fopen('population_data.dat', 'w');
fwrite(fid, 'y0');
fclose(fid);
fid = fopen('population_data.dat','r');
datacell = textscan(fid, '%f%f%f', 'HeaderLines', 1, 'Collect', 1);
csvwrite('population_data.csv',y0)
type 'population_data.dat'

%X_matrix = (y0 - min(y0,[],2))./(max(y0,[],2) - min(y0,[],2));
