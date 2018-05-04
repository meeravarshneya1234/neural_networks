settings.model_name = 'Grandi';
%modelnames = {'TT06','Ohara'};
% options - 'Fox', 'Hund', 'Livshitz',
% 'Devenyi','Shannon','TT04','TT06','Grandi','Ohara'
% 'TT06_opt','Grandi_opt','Ohara_opt'

settings.celltype = 'endo'; % size should be same as model_name, enter one for each model
% options only available for human models - 'epi', 'endo', 'mid',

settings.PCL =1000 ;  % Interval bewteen stimuli,[ms]
settings.stim_delay = 100 ; % Time the first stimulus, [ms]
settings.stim_dur = 2 ; % Stimulus duration
settings.stim_amp = 20.6; % Stimulus amplitude 
settings.nBeats = 100 ; % Number of beats to simulate 
settings.numbertokeep =1;% Determine how many beats to keep. 1 = last beat, 2 = last two beats
settings.steady_state = 1;

% Ks and Kr scaling must have the same length vector. Mainly change
% Ks_scale and Kr_scale when you want to test what different ratios of the
% two look like 
settings.Ks_scale = 1; % perturb IKs, set to 1 for baseline
settings.Kr_scale = 1; % perturb IKr, set to 1 for baseline
Ca_scale = [1 1.5]; % perturb ICaL, set to 1 for baseline

settings.variations = 5000;
settings.sigma = 0.2;
settings.scalings = exp(settings.sigma*randn(18,settings.variations))' ; % same parameter variation for each pop

%%     
str = {'ctl','exp'};
for i = 1:length(Ca_scale) 
    settings.Ca_scale = Ca_scale(i);
    X = pop_program(settings);  
    datatable.(str{i}) = X;  
end 

%%
for i = 1:length((str))
    figure
    for ii = 1:settings.variations % only plot the first 20 for the figure
        plot(datatable.(str{i}).times{ii},datatable.(str{i}).V{ii},'linewidth',2)
        hold on
    end
    title(str{i})
    xlabel('time (ms)')
    ylabel('voltage (mV)')
    ylim([-100 60])
      
end

%%
modelnames = {'ctl','exp'};
t_cutoffs = [4 4];
flags =[0 0 ];
for i = 1:length(modelnames)
    %load([modelnames{i} 'pop.mat'])
    
    x = datatable.(str{i}).times(:,1);
    y = datatable.(str{i}).V(:,1);
    
    figure
    subplot(3,1,1)
    hold on
    cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
    title(modelnames{i})
    
    flag = flags(i);
    t_cutoff = t_cutoffs(i);
    [APfails,EADs] = cleandata(datatable.(str{i}).APDs(:,1),datatable.(str{i}).times(:,1),datatable.(str{i}).V(:,1),t_cutoff,flag);
    [number_of_failed,ind_failed] = find(APfails ==1); % number of failed to repolarize
    [number_of_EADs,ind_EADs] = find(EADs==1); % number of EADs
    
    indexs = [ind_EADs ind_failed];
    x = datatable.(str{i}).times(indexs,1);
    y = datatable.(str{i}).V(indexs,1);
    
    figure(gcf)
    subplot(3,1,2)
    hold on
    cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
    title(modelnames{i})
    
    clean_datatable = [];
    
    clean_datatable.times = datatable.(str{i}).times(~(EADs' + APfails'));
    clean_datatable.V= datatable.(str{i}).V(~(EADs' + APfails'));
    clean_datatable.states = datatable.(str{i}).states(~(EADs' + APfails'));
    clean_datatable.APDs = datatable.(str{i}).APDs(~(EADs' + APfails'));
    clean_datatable.scaling = datatable.(str{i}).scalings(~(EADs' + APfails'),:);
    
    x = clean_datatable.times(:,1);
    y = clean_datatable.V(:,1);
    
    figure(gcf)
    subplot(3,1,3)
    hold on
    cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
    title([modelnames{i} ' noEADs'])  
    withEADs(i) = 5000 - length(x);
%     
%     if noEADs(i) > 0 
%        settings.model_name = modelnames{i};
%        settings.celltype = celltypes{i};
%        settings.stim_amp = amps(i);
%        settings.variations = noEADs(i);
%        X = rerunAPs(settings,clean_datatable);
%        save([modelnames{i} 'cleanpop.mat'],'X')
%        
%        x = X.times(:,1);
%        y = X.V(:,1);
%        
%        figure
%        hold on
%        cellfun(@(x,y) plot(x,y,'linewidth',2),x,y)
%        title([modelnames{i} 'cleaned'])      
%     end 
    
end 
%%
% for i = 1:length(clean_datatable.times(:,1))
%     figure
%     h = axes;
%     plot(clean_datatable.times{i,1},clean_datatable.V{i,1},'k','linewidth',2);
%     %ylim([-100 50])
%     set(gcf,'units','pixels','position',[489,398,281 253])
%     set(h,'Visible','off')
%     saveas(gcf,[num2str(i) '.jpg'])    
%     close(gcf)
% end 
