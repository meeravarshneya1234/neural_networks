% Updated Vrest calculation method and inclusion of drugs that caused EADs



%% Load V, t, and Ca traces for each simulation

path = '/Users/Jake/Dropbox/jake sobie lab/MATLAB code for project/Megan_sim/';
%% load drug groups and outputs from all conditions

load([path,'drugdata_k2013_m2011'],'datacell');
names = datacell(2:end,1);
group = cell2mat(datacell(2:end,2)); %get categories from datacell
clear datacell

n_drugs = length(names);

%% Locations of all simulation data
%9 conditions: endo, mid, and epi each at 2000, 1000, and 500 ms pacing intervals
paths = {'SimulationData/ORd_01xD/ORd.epi.01xD.2000_0/Drugs_ORd_01xD_2000_0',[];...  %endo 2000
         'SimulationData/ORd_01xD/ORd.epi.01xD.1000_0/Drugs_ORd_01xD_1000_0',[];...  %endo 1000
         'SimulationData/ORd_01xD/ORd.epi.01xD.500_0/Drugs_ORd_01xD_500_0',[];...    %endo 500
         'SimulationData/ORd_01xD/ORd.epi.01xD.2000_2/Drugs_ORd_01xD_2000_2',[];...   %mid 2000
         'SimulationData/ORd_01xD/ORd.epi.01xD.1000_2/Drugs_ORd_01xD_1000_2',[];...   %mid 1000
         'SimulationData/ORd_01xD/ORd.epi.01xD.500_2/Drugs_ORd_01xD_500_2',[];...     %mid 500
         'SimulationData/ORd_01xD/ORd.epi.01xD.2000_1/Drugs_ORd_01xD_2000_1',[];...   %epi 2000
         'SimulationData/ORd_01xD/ORd.epi.01xD.1000_1/Drugs_ORd_01xD_1000_1',[];...   %epi 1000
         'SimulationData/ORd_01xD/ORd.epi.01xD.500_1/Drugs_ORd_01xD_500_1',[]};       %epi 500

%% Load baseline measurements needed for later calcs
load([path,'SimulationData/outputs_bl'],'outputs_bl')


%% Calculate outputs under each condition
outputs_all = [];

% load replace %replace is a 9x86 matrix (conditions x drugs) which indicates which drugs results need to be replaced under which conditions
% load paths_replace

for i = 1:size(paths,1)  %loop through conditions
    load([path,paths{i,1}],'actionpotentials','times', 'catransients', 'finalbeat_ind')
    aps = actionpotentials;
    cats = catransients;
    T = times;
    fbi = finalbeat_ind;
    if paths{i,2}       %if two drugs sets were simulated separately combine them
       load([path,paths{i,2}],'actionpotentials','times', 'catransients', 'finalbeat_ind')
       aps = [aps;actionpotentials];
       cats = [cats;catransients];
       T = [T;times];
       fbi = [fbi;finalbeat_ind];
    end
    
    
    %replace drugs that did not repolarize/had EADs at EFTPC with
    %measurements from lower doses 
%     if any(replace(i,:))
%         drugs = find(replace(i,:));
%         for iii= drugs
%             load([path,paths_replace{i,iii}],'actionpotentials','times', 'catransients', 'finalbeat_ind')
%             if length(actionpotentials)>1
%                 aps(iii) = actionpotentials(iii);
%                 cats(iii) = catransients(iii);
%                 T(iii) = times(iii);
%                 fbi(iii) = finalbeat_ind(iii);
%             else    
%                 aps(iii) = actionpotentials;
%                 cats(iii) = catransients;
%                 T(iii) = times;
%                 fbi(iii) = finalbeat_ind;
%             end
%         end    
%     end
    %% Calculate outputs
    
    outputnames = {'Resting V';'Max dVdt';'Peak V';'APD';'APD50';'APD90';...
    'AP triangulation';'CaT amp';'Peak Ca';'CaD50';'CaD90';'CaT triangulation';'dCa'};

    Vrest = zeros(1,n_drugs);
    upstrokes = zeros(1,n_drugs);
    Vpeak = zeros(1,n_drugs);
    APD = zeros(1,n_drugs);
    APD50 = zeros(1,n_drugs);
    APD90 = zeros(1,n_drugs);
    TriAP = zeros(1,n_drugs);
    DCai = zeros(1,n_drugs);
    Capeak = zeros(1,n_drugs);
    CaD50 = zeros(1,n_drugs);
    CaD90 = zeros(1,n_drugs);
    TriCa = zeros(1,n_drugs);
    dCa = zeros(1,n_drugs);
    
    for ii = 1:n_drugs    %loop through drugs
        beatstart = fbi{ii};
        V = aps{ii}(beatstart:end);
        Cai = cats{ii}(beatstart:end);
        t = T{ii}(beatstart:end);

        %% Calculate outputs
        % Vrest, Vpeak, APD, Vrest, dCai, APD90, APD50, triangulation (90-50), max upstroke,
        % Ca duration at 90 and at 50, Ca triangulation, Peak Ca

        %Voltage outputs   
        Vderiv = diff(V)./diff(t) ;
        [dVdtmax,dexmax] = max(Vderiv) ;
        tinit = t(dexmax) ; %Time of maximum dV/dt, consider this beginning of action potential

%         Vrest(ii) = V(1) ;
        vrest = min(V(1:dexmax));
 
        %%Then determine peak voltage of action potential, for two reasons
        %%1) Because repolarization must, by definition, occur after this
        %%2) To compute 50%, 90%, etc., must have this value
        [peakV,peakdex] = max(V) ;        
        tpeak = t(peakdex) ;
        repoldex = find(t > tpeak & V < -60) ;
        V50_exact = 0.5*(peakV - vrest) + vrest ;
        V90_exact = 0.1*(peakV - vrest) + vrest ;
        V50dex = find(t > tpeak & V < V50_exact) ;
        V90dex = find(t > tpeak & V < V90_exact);
        
        DCai(ii) = max(Cai)-min(Cai);
        [peakCa,peakdex] = max(Cai) ;
        Capeak(ii) = peakCa ;
        tpeak_Ca = t(peakdex) ;
        Ca90_exact = 0.1*(DCai(ii)) + min(Cai) ;
        Ca90dex = find(t > tpeak_Ca & Cai < Ca90_exact);
        
        if any(V90dex) && any(repoldex) && vrest < -60 && any(Ca90dex) %trace is good to use: repolarized before and after
            Vpeak(ii) = peakV ;
            Vrest(ii) = vrest;
            upstrokes(ii) = dVdtmax ;
            repoltime = t(repoldex(1)) ;
            APD(ii) = repoltime - tinit ;
            V50_time = t(V50dex(1));
            APD50(ii) = V50_time - tinit;
            V90_time = t(V90dex(1));
            APD90(ii) = V90_time - tinit;
            TriAP(ii) = V90_time-V50_time; %Triangulation

            %% Ca outputs
%             DCai(ii) = max(Cai)-min(Cai);
%             [peakCa,peakdex] = max(Cai) ;
%             Capeak(ii) = peakCa ;
%             tpeak = t(peakdex) ;
            %CaT duration calcs
            Ca50_exact = 0.5*(DCai(ii)) + min(Cai) ; 
%             Ca90_exact = 0.1*(DCai(ii)) + min(Cai) ;
            Ca50dex = find(t > tpeak_Ca & Cai < Ca50_exact) ;
            Ca50_time = t(Ca50dex(1));
            CaD50(ii) = Ca50_time - tinit;
%             Ca90dex = find(t > tpeak & Cai < Ca90_exact);
            Ca90_time = t(Ca90dex(1));
            CaD90(ii) = Ca90_time - tinit;
            TriCa(ii) = Ca90_time-Ca50_time;                 %Triangulation

            dCa(ii) = Cai(1);
            
        else %check if 2-to-last interval repolarizes before last stim
            V = aps{ii}(1:beatstart+dexmax);
            Cai = cats{ii}(1:beatstart+dexmax);
            t = T{ii}(1:beatstart+dexmax);

            
            %% Calculate outputs
            % Vrest, Vpeak, APD, Vrest, dCai, APD90, APD50, triangulation (90-50), max upstroke,
            % Ca duration at 90 and at 50, Ca triangulation, Peak Ca

            %Voltage outputs   
            Vderiv = diff(V)./diff(t) ;
            [dVdtmax,dexmax] = max(Vderiv(1:beatstart)) ;  %interval includes small amount of next upstroke, don't want to use this so end at "beatstart"
            tinit = t(dexmax) ; %Time of maximum dV/dt, consider this beginning of action potential

    %         Vrest(ii) = V(1) ;
            vrest = min(V(1:dexmax));

            %%Then determine peak voltage of action potential, for two reasons
            %%1) Because repolarization must, by definition, occur after this
            %%2) To compute 50%, 90%, etc., must have this value
            [peakV,peakdex] = max(V) ;        
            tpeak = t(peakdex) ;
            repoldex = find(t > tpeak & V < -60) ;
            V50_exact = 0.5*(peakV - vrest) + vrest ;
            V90_exact = 0.1*(peakV - vrest) + vrest ;
            V50dex = find(t > tpeak & V < V50_exact) ;
            V90dex = find(t > tpeak & V < V90_exact);
            
            DCai(ii) = max(Cai)-min(Cai);
            [peakCa,peakdex] = max(Cai) ;
            Capeak(ii) = peakCa ;
            tpeak_Ca = t(peakdex) ;
            Ca90_exact = 0.1*(DCai(ii)) + min(Cai) ;
            Ca90dex = find(t > tpeak_Ca & Cai < Ca90_exact);

            
            if any(V90dex) && any(repoldex) && vrest < -60 && any(Ca90dex) %trace is good to use: repolarized before and after
                disp(['Using 2nd-to-last for ',num2str(ii),'_',num2str(i)])
                Vpeak(ii) = peakV ;
                Vrest(ii) = vrest;
                upstrokes(ii) = dVdtmax ;
                repoltime = t(repoldex(1)) ;
                APD(ii) = repoltime - tinit ;
                V50_time = t(V50dex(1));
                APD50(ii) = V50_time - tinit;
                V90_time = t(V90dex(1));
                APD90(ii) = V90_time - tinit;
                TriAP(ii) = V90_time-V50_time; %Triangulation

                %% Ca outputs
%                 DCai(ii) = max(Cai)-min(Cai);
%                 [peakCa,peakdex] = max(Cai) ;
%                 Capeak(ii) = peakCa ;
%                 tpeak = t(peakdex) ;
                %CaT duration calcs
                Ca50_exact = 0.5*(DCai(ii)) + min(Cai) ; 
%                 Ca90_exact = 0.1*(DCai(ii)) + min(Cai) ;
                Ca50dex = find(t > tpeak_Ca & Cai < Ca50_exact) ;
                Ca50_time = t(Ca50dex(1));
                CaD50(ii) = Ca50_time - tinit;
%                 Ca90dex = find(t > tpeak & Cai < Ca90_exact);
                Ca90_time = t(Ca90dex(1));
                CaD90(ii) = Ca90_time - tinit;
                TriCa(ii) = Ca90_time-Ca50_time;                 %Triangulation
 
                dCa(ii) = Cai(1);
            else   %dose is no good 
                disp([num2str(ii),'_',num2str(i),' does not repolarize'])
                Vpeak(ii) = NaN ;
                Vrest(ii) = NaN;
                upstrokes(ii) = NaN;
                APD(ii) = NaN;
                APD50(ii) = NaN;
                APD90(ii) = NaN;
                TriAP(ii) = NaN;
                CaD50(ii) = NaN;
                CaD90(ii) = NaN;
                TriCa(ii) = NaN;
                DCai(ii) = NaN;
                Capeak(ii) = NaN ;
                dCa(ii) = NaN;
            end
        end
    end

    outputs = [Vrest;upstrokes;Vpeak;APD;APD50;APD90;TriAP;DCai;Capeak;...  %12x86
        CaD50;CaD90;TriCa;dCa];
    outputs_all = [outputs_all;outputs];
    
end

%Calculate TDR, RD1, RD2
diff(1,:) = outputs_all(17,:) - outputs_all(95,:);   %APD endo1000 - epi1000
diff(2,:) = outputs_all(17,:) - outputs_all(56,:);   %APD endo1000 - mid1000
diff(3,:) = outputs_all(95,:) - outputs_all(56,:);   %APD epi1000  - mid1000 

TDR = max(abs(diff),[],1);

RD1_endo = outputs_all(4,:)-outputs_all(30,:);   %APD endo2000 - endo500
RD2_endo = abs((outputs_all(4,:)-outputs_bl(4))./outputs_bl(4))-...
    abs((outputs_all(30,:)-outputs_bl(30))./outputs_bl(30));

RD1_mid = outputs_all(43,:)-outputs_all(69,:); %APD mid2000 - mid500
RD2_mid = abs((outputs_all(43,:)-outputs_bl(43))./outputs_bl(43))-...
    abs((outputs_all(69,:)-outputs_bl(69))./outputs_bl(69));

RD1_epi = outputs_all(82,:)-outputs_all(108,:); %APD epi2000 - epi500
RD2_epi = abs((outputs_all(82,:)-outputs_bl(82))./outputs_bl(82))-...
    abs((outputs_all(108,:)-outputs_bl(108))./outputs_bl(108));

outputs_all = [outputs_all;RD1_endo;RD1_mid;RD1_epi;RD2_endo;RD2_mid;RD2_epi;TDR]';  


%% Output names for each pacing rate and cell type
freq = {' 2000',' 1000',' 500'};
outputnames_endo = {};
outputnames_mid = {};
outputnames_epi = {};
for ii = 1:length(freq)
    for i = 1:length(outputnames)
       nextname_endo = [outputnames{i},freq{ii},' endo']; 
       outputnames_endo = [outputnames_endo,nextname_endo]; %1x36 cell
       nextname_mid = [outputnames{i},freq{ii},' mid']; 
       outputnames_mid = [outputnames_mid,nextname_mid]; %1x36 cell
       nextname_epi = [outputnames{i},freq{ii},' epi']; 
       outputnames_epi = [outputnames_epi,nextname_epi]; %1x36 cell
    end
end
outputnames_all = [outputnames_endo,outputnames_mid,outputnames_epi,...
    'RD1_endo','RD1_mid','RD1_epi','RD2_endo','RD2_mid','RD2_epi','TDR',]; %1x115 cell


%% Remove drugs with NaN outputs
[remove,~] = find(isnan(outputs_all));
remove = unique(remove);
outputs_all(remove,:) = [];
names(remove) = []; 
group(remove) = [];
n_drugs = n_drugs-length(remove);


%% Plot heatmap of outputs for all drugs
%don't include RD2 because already a percent change
% 121:123
outputs_fold = log(abs([outputs_all(:,1:120),outputs_all(:,124)]...
    ./(ones(n_drugs,1)*...
    [outputs_bl(:,1:120),outputs_bl(:,124)])))...
    /log(2); %log base 2 of ratio

%clustergram function clusters both axes hierarchically and generates figure with heatmap and dendrogram
clustergram(outputs_fold','RowLabels',[outputnames_all(1:120),outputnames_all(124)],...
    'ColumnLabels',names,'Colormap','jet')


%% Output correlations
corrmat = corrcoef(outputs_all);

%order outputs based on hierarchical clustering
clustergram(corrmat,'RowLabels',outputnames_all,'ColumnLabels',outputnames_all,...
'Colormap','jet')

save 'outputs' 'outputs'
save 'outputnames' 'outputnames'