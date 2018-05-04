n = 


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
    