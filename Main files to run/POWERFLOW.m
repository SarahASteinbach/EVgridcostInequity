%% Run power flow and determine overloads



%% Power FLOW RURAL CASE
% case selector:
% 1: charging loads induced by homecharging EVs
% 2: household electrical load profile
% 3: combined profiles of charging loads and households


%Adjust EV penetration based on income group analyzed
%EV_pen=22.4;
EV_pen=35.7;
%EV_pen=31.1;


iterations = 800;
timesteps = 2016;         
time_vars = [december,march,june,september];
% Only explore december as most challenging season
number_of_seasons = 1;
rural.no_nodes = 95;


%Technical model specifications (as code of a more flexible model used)
selected_cases = [3];
house_scaling_factor=0;


parfor time=1:number_of_seasons
    time_of_year = time_vars(time);
    [initial_EV_socs] = simulateSOCs(rural,time_of_year,housetype_charger_available, H16_car,car_segment_distribution, car_data, samples, speed,timestep_size, no_timesteps);
for EV_penetration = EV_pen
    min_voltage = zeros(timesteps,iterations*length(selected_cases));
   
    min_voltageID = zeros(timesteps,iterations*length(selected_cases));
    overloaded_branchID = zeros(timesteps,iterations*length(selected_cases));

    total_load = zeros(timesteps,iterations*length(selected_cases));
    total_HH_load = zeros(timesteps,iterations*length(selected_cases));
    total_EV_load = zeros(timesteps,iterations*length(selected_cases));
    branch_log = cell(timesteps,iterations*length(selected_cases));
    CriticalUsefactor=zeros(timesteps,iterations*length(selected_cases));
    CriticalBranchID=zeros(timesteps,iterations*length(selected_cases));
    TransformerUsefactor= zeros(timesteps,iterations*length(selected_cases));
    crit_fact_plot = nan(timesteps,iterations*length(selected_cases));
    hh_load_plot = nan(timesteps,iterations*11);
    hh_load_plot_idx = 1;
    for k = 1:iterations
        [nodes, charging_nodes] = assembleNodesN1(rural,time_of_year,EV_penetration,housetype_charger_available, H16_car,car_segment_distribution, car_data, samples, speed,timestep_size, no_timesteps,initial_EV_socs,house_scaling_factor);
        for l=1:length(selected_cases)
            selected_case=selected_cases(l);
            mpc = loadcase('SimBench_Rural2');
            %run load flow analysis without printing anything and enfore reactive limits
            mpopt = mpoption('pf.nr.max_it',100,'verbose',0,'out.all',0, 'pf.enforce_q_lims', 1);

            for i=1:timesteps
                % skip transformer busses für rural2 grid
                % skip bus 104 (n=97) and bus 19
                for n = 1:95
                    %skip bus 19 
                    if n < 19
                        % get profiles and convert in MW
                        profile = nodes{n,selected_case}/1000;
                        total_HH_load(i,iterations*(l-1)+k) = total_HH_load(i,iterations*(l-1)+k) + nodes{n,2}(i)/1000;
                        total_EV_load(i,iterations*(l-1)+k) = total_EV_load(i,iterations*(l-1)+k) + nodes{n,1}(i)/1000;
                        % Assign power demand to bus
                            mpc.bus(n,2) = 1;
                            mpc.bus(n,3) = -profile(i)*0.93;
                            mpc.bus(n,4) = sqrt((profile(i))^2-(profile(i)*0.93)^2);
                            % set gen                
                            mpc.gen(n,2) = 0;
                            mpc.gen(n,9) = 0;                
                    else
                        % get profiles and convert in MW
                        profile = nodes{n,selected_case}/1000;
                        total_HH_load(i,iterations*(l-1)+k) = total_HH_load(i,iterations*(l-1)+k) + nodes{n,2}(i)/1000;
                        total_EV_load(i,iterations*(l-1)+k) = total_EV_load(i,iterations*(l-1)+k) + nodes{n,1}(i)/1000;

                        % if negative act as power demand
                        % Assign power demand to bus
                            mpc.bus(n+1,2) = 1;
                            mpc.bus(n+1,3) = -profile(i)*0.93;
                            mpc.bus(n+1,4) = sqrt((profile(i))^2-(profile(i)*0.93)^2);
                            % set gen                
                            mpc.gen(n+1,2) = 0;
                            mpc.gen(n+1,9) = 0;                
                    end
                end   
                result = runpf(mpc,mpopt);
                %Find critical branches
                Abranch=abs(complex(result.branch(2:end,14),result.branch(2:end,15))); %creates vector with apparent power values at the branches
                Capacity=mpc.branch(2:end,6); %creates a vector with the capity values of the branches
                Usefactor=Abranch./Capacity; %creates a vector with the usefactors at the branches
                CriticalUsefactor(i,iterations*(l-1)+k)=max(Usefactor); %determines the critical usefactor
                CriticalBranchID(i,iterations*(l-1)+k)=find(Usefactor == max(Usefactor(:)),1,'first'); %determines the critical branch
               
                if (result.branch(1,14)<0)
                    TransformerUsefactor(i,iterations*(l-1)+k) = result.branch(1,14)./mpc.branch(1,6);
                else
                    TransformerUsefactor(i,iterations*(l-1)+k) = abs(complex(result.branch(1,14),result.branch(1,15)))./mpc.branch(1,6);
                end
                %Min voltage violations were analyzed but did not occur,
                %hence not mentioned in final paper
                min_voltage(i,iterations*(l-1)+k) = min(result.bus(:, 8));
                
                min_voltageID(i,iterations*(l-1)+k) = find(result.bus(:, 8) == min(result.bus(:, 8)),1,'first');
                %only if usefactor >1 or min_voltage<0.9
                if CriticalUsefactor(i,iterations*(l-1)+k) > 1
                    overloaded_branchID(i,iterations*(l-1)+k) = CriticalBranchID(i,iterations*(l-1)+k);
                elseif min_voltageID(i,iterations*(l-1)+k) < 0.9
                    overloaded_branchID(i,iterations*(l-1)+k) = min_voltageID(i,iterations*(l-1)+k);
                elseif TransformerUsefactor(i,iterations*(l-1)+k) > 1
                    overloaded_branchID(i,iterations*(l-1)+k) = 1;                             
                end
                
                total_load(i,iterations*(l-1)+k) = sum(result.bus(:, 3));
                branch_log{i,1} = result.branch(:,:);   
            end
        end
        disp(k)
    end
    for c=1:length(selected_cases)
        selected_case=selected_cases(c);
        data_plot = CriticalUsefactor(1:timesteps,iterations*(c-1)+1:iterations*c);
        h1 = figure('visible','off');
        hold on
        plot(ones(timesteps))
        fanChart(1:size(data_plot,1), data_plot);
        ylim([0 1.5])
        xlim([0 timesteps])
        xlabel('Time') 
        ylabel('Load level (in %)')
        set(gca,'xtick',0:144:timesteps)
        set(gca,'xticklabel',{'Mon 0:00' 'Mon 12:00' 'Tue 0:00' 'Tue 12:00' 'Wed 0:00' 'Wed 12:00' 'Thu 0:00' 'Thu 12:00' 'Fri 0:00' 'Fri 12:00' 'Sat 0:00' 'Sat 12:00' 'Sun 0:00' 'Sun 12:00' })
        set(gca,'fontname','times', 'FontSize', 12)
        set(gcf,'position',[0,0,1200,1000])

        %legend('Pct10','Pct20','Pct30','Pct40','Pct50','Pct60','Pct70','Pct80','Pct90','Median','Location','southoutside')
        saveas(h1,sprintf('Results/Rural/%s Crit Branch Loads Case %d EV_penetration %d.png',time_of_year.name,selected_cases(c),EV_penetration));
        close(h1);
        hold off
        %save critical branch usefactor as csv
        writematrix(CriticalUsefactor(1:timesteps,iterations*(c-1)+1:iterations*c),sprintf('Results/Rural/%s Critical Branch Loads %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));


        h2=figure('visible','off');
        hold on
        plot(-mean(total_HH_load(:,iterations*(c-1)+1:iterations*c),2))
        plot(-mean(total_HH_load(:,iterations*(c-1)+1:iterations*c),2)-mean(total_EV_load(:,iterations*(c-1)+1:iterations*c),2))
        xlim([0 timesteps])
        xlabel('Time') 
        ylabel('Load (in MWh)')
        set(gca,'xtick',0:144:timesteps)
        set(gca,'xticklabel',{'Mon 0:00' 'Mon 12:00' 'Tue 0:00' 'Tue 12:00' 'Wed 0:00' 'Wed 12:00' 'Thu 0:00' 'Thu 12:00' 'Fri 0:00' 'Fri 12:00' 'Sat 0:00' 'Sat 12:00' 'Sun 0:00' 'Sun 12:00' })
        set(gca,'fontname','times', 'FontSize', 12)
        legend('Total Household Loads','Total Household and EV Loads', 'Location','southoutside')
        
        hold off
        saveas(h2,sprintf('Results/Rural/%s Loads Case %d EV Penetration %d.png',time_of_year.name, selected_case,EV_penetration));
        close(h2);
        if selected_case == 3 
            %save HH, EV as csv
            current_loads=[total_HH_load(:,iterations*(c-1)+1:iterations*c), total_HH_load(:,iterations*(c-1)+1:iterations*c) + total_EV_load(:,iterations*(c-1)+1:iterations*c)];
            writematrix(current_loads,sprintf('Results/Rural/%s HH Loads %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));
        end



        %plot transformer load
        h3=figure('visible','off');
        hold on
        data_plot = TransformerUsefactor(1:timesteps,iterations*(c-1)+1:iterations*c);
        plot(ones(timesteps))
        fanChart(1:size(data_plot,1), data_plot);
        ylim([0 1.5])
        xlim([0 timesteps])
        xlabel('Time') 
        ylabel('Load level (in %)')
        set(gca,'xtick',0:144:timesteps)
        set(gca,'xticklabel',{'Mon 0:00' 'Mon 12:00' 'Tue 0:00' 'Tue 12:00' 'Wed 0:00' 'Wed 12:00' 'Thu 0:00' 'Thu 12:00' 'Fri 0:00' 'Fri 12:00' 'Sat 0:00' 'Sat 12:00' 'Sun 0:00' 'Sun 12:00' })
        set(gca,'fontname','times', 'FontSize', 12)
        set(gcf,'position',[0,0,1200,1000])
        %legend('Pct10','Pct20','Pct30','Pct40','Pct50','Pct60','Pct70','Pct80','Pct90','Median','Location','southoutside')
        saveas(h3,sprintf('Results/Rural/%s Transformer Loads Case %d EV_penetration %d.png',time_of_year.name,selected_cases(c),EV_penetration));
        close(h3);
        %save critical transformer usefactor as csv
        writematrix(TransformerUsefactor(1:timesteps,iterations*(c-1)+1:iterations*c),sprintf('Results/Rural/%s Critical Transformer Loads %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));


        %plot min branch voltage
        h4=figure('visible','off');
        hold on
        data_plot = min_voltage(1:timesteps,iterations*(c-1)+1:iterations*c);
        plot(0.9*ones(timesteps))
        fanChart(1:size(data_plot,1), data_plot);
        ylim([0.8 1.05])
        xlim([0 timesteps])
        xlabel('Time') 
        ylabel('Voltage level (in %)')
        set(gca,'xtick',0:144:timesteps)
        set(gca,'xticklabel',{'Mon 0:00' 'Mon 12:00' 'Tue 0:00' 'Tue 12:00' 'Wed 0:00' 'Wed 12:00' 'Thu 0:00' 'Thu 12:00' 'Fri 0:00' 'Fri 12:00' 'Sat 0:00' 'Sat 12:00' 'Sun 0:00' 'Sun 12:00' })
        set(gca,'fontname','times', 'FontSize', 12)
        set(gcf,'position',[0,0,1200,1000])
        %legend('Pct10','Pct20','Pct30','Pct40','Pct50','Pct60','Pct70','Pct80','Pct90','Median','Location','southoutside')
        saveas(h4,sprintf('Results/Rural/%s Min voltage Case %d EV_penetration %d.png',time_of_year.name,selected_cases(c),EV_penetration));
        close(h4);
        close all;
        %save min voltage as csv
        writematrix(min_voltage(1:timesteps,iterations*(c-1)+1:iterations*c),sprintf('Results/Rural/%s Min voltage %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));


        %save critical branch IDs as csv
        writematrix(overloaded_branchID(:,iterations*(c-1)+1:iterations*c),sprintf('Results/Rural/%s Critical Branches %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));
    end
end
end





% Power FLOW SEMIURBAN CASE
% case selector:
% 1: charging loads induced by homecharging EVs
% 2: household electrical load profile
% 3: combined profiles of charging loads and households

semi_urb.no_nodes = 109;

parfor time=1:number_of_seasons
    time_of_year = time_vars(time);
    [initial_EV_socs] = simulateSOCs(semi_urb,time_of_year,housetype_charger_available, H16_car,car_segment_distribution, car_data, samples, speed,timestep_size, no_timesteps);
for EV_penetration = EV_pen
    min_voltage = zeros(timesteps,iterations*length(selected_cases));
    
    min_voltageID = zeros(timesteps,iterations*length(selected_cases));
    overloaded_branchID = zeros(timesteps,iterations*length(selected_cases));

    total_load = zeros(timesteps,iterations*length(selected_cases));
    total_HH_load = zeros(timesteps,iterations*length(selected_cases));
    total_EV_load = zeros(timesteps,iterations*length(selected_cases));
    branch_log = cell(timesteps,iterations*length(selected_cases));
    CriticalUsefactor=zeros(timesteps,iterations*length(selected_cases));
    CriticalBranchID=zeros(timesteps,iterations*length(selected_cases));
    TransformerUsefactor= zeros(timesteps,iterations*length(selected_cases));
    crit_fact_plot = nan(timesteps,iterations*length(selected_cases));
    hh_load_plot = nan(timesteps,iterations*11);
    hh_load_plot_idx = 1;
    for k = 1:iterations
        [nodes, charging_nodes] = assembleNodesN1(semi_urb,time_of_year,EV_penetration,housetype_charger_available, H16_car,car_segment_distribution, car_data, samples, speed,timestep_size, no_timesteps,initial_EV_socs,house_scaling_factor);
        for l=1:length(selected_cases)
            selected_case=selected_cases(l);
            mpc = loadcase('SimBench_SemiUrb5');
            %run load flow analysis without printing anything and enfore reactive limits
            mpopt = mpoption('pf.nr.max_it',100,'verbose',0,'out.all',0, 'pf.enforce_q_lims', 1);

            for i=1:timesteps
                % skip transformer busses of SemiUrb5 Grid
                % skip bus 114 (n = 111) and skip bus 73 (n = 71)
                for n = 1:109
                    % skip bus 73 (n = 71)
                    if n < 71
                        % get profiles and convert in MW
                        profile = nodes{n,selected_case}/1000;
                        total_HH_load(i,iterations*(l-1)+k) = total_HH_load(i,iterations*(l-1)+k) + nodes{n,2}(i)/1000;
                        total_EV_load(i,iterations*(l-1)+k) = total_EV_load(i,iterations*(l-1)+k) + nodes{n,1}(i)/1000;
                        % if negative act as power demand
                             % set bus
                            mpc.bus(n,2) = 1;
                            mpc.bus(n,3) = -profile(i)*0.93;
                            mpc.bus(n,4) = sqrt((profile(i))^2-(profile(i)*0.93)^2);
                            % set gen                
                            mpc.gen(n,2) = 0;
                            mpc.gen(n,9) = 0;                
                         
                    else
                        % get profiles and convert in MW
                        profile = nodes{n,selected_case}/1000;
                        total_HH_load(i,iterations*(l-1)+k) = total_HH_load(i,iterations*(l-1)+k) + nodes{n,2}(i)/1000;
                        total_EV_load(i,iterations*(l-1)+k) = total_EV_load(i,iterations*(l-1)+k) + nodes{n,1}(i)/1000;
                        
                            % set bus
                            mpc.bus(n+1,2) = 1;
                            mpc.bus(n+1,3) = -profile(i)*0.93;
                            mpc.bus(n+1,4) = sqrt((profile(i))^2-(profile(i)*0.93)^2);
                            % set gen                
                            mpc.gen(n+1,2) = 0;
                            mpc.gen(n+1,9) = 0;                
                        end
                end   
                result = runpf(mpc,mpopt);
                %Find critical branches
                Abranch=abs(complex(result.branch(2:end,14),result.branch(2:end,15))); %creates vector with apparent power values at the branches
                Capacity=mpc.branch(2:end,6); %creates a vector with the capity values of the branches
                Usefactor=Abranch./Capacity; %creates a vector with the usefactors at the branches
                CriticalUsefactor(i,iterations*(l-1)+k)=max(Usefactor); %determines the critical usefactor
                CriticalBranchID(i,iterations*(l-1)+k)=find(Usefactor == max(Usefactor(:)),1,'first'); %determines the critical branch
               
                if (result.branch(1,14)<0)
                    TransformerUsefactor(i,iterations*(l-1)+k) = result.branch(1,14)./mpc.branch(1,6);
                else
                    TransformerUsefactor(i,iterations*(l-1)+k) = abs(complex(result.branch(1,14),result.branch(1,15)))./mpc.branch(1,6);
                end
                %Min voltage violations were analyzed but did not occur,
                %hence not mentioned in final paper
                min_voltage(i,iterations*(l-1)+k) = min(result.bus(:, 8));
                
                min_voltageID(i,iterations*(l-1)+k) = find(result.bus(:, 8) == min(result.bus(:, 8)),1,'first');
                %only if usefactor >1 or min_voltage<0.9
                if CriticalUsefactor(i,iterations*(l-1)+k) > 1
                    overloaded_branchID(i,iterations*(l-1)+k) = CriticalBranchID(i,iterations*(l-1)+k);
                elseif min_voltageID(i,iterations*(l-1)+k) < 0.9
                    overloaded_branchID(i,iterations*(l-1)+k) = min_voltageID(i,iterations*(l-1)+k);
                elseif TransformerUsefactor(i,iterations*(l-1)+k) > 1
                    overloaded_branchID(i,iterations*(l-1)+k) = 1;             
                end
                
                total_load(i,iterations*(l-1)+k) = sum(result.bus(:, 3));
                branch_log{i,1} = result.branch(:,:);   
            end
        end
        disp(k)
    end
    for c=1:length(selected_cases)
        selected_case=selected_cases(c);
        data_plot = CriticalUsefactor(1:timesteps,iterations*(c-1)+1:iterations*c);
        h1 = figure('visible','off');
        hold on
        plot(ones(timesteps))
        fanChart(1:size(data_plot,1), data_plot);
        ylim([0 1.5])
        xlim([0 timesteps])
        xlabel('Time') 
        ylabel('Load level (in %)')
        set(gca,'xtick',0:144:timesteps)
        set(gca,'xticklabel',{'Mon 0:00' 'Mon 12:00' 'Tue 0:00' 'Tue 12:00' 'Wed 0:00' 'Wed 12:00' 'Thu 0:00' 'Thu 12:00' 'Fri 0:00' 'Fri 12:00' 'Sat 0:00' 'Sat 12:00' 'Sun 0:00' 'Sun 12:00' })
        set(gca,'fontname','times', 'FontSize', 12)
        set(gcf,'position',[0,0,1200,1000])
        %legend('Pct10','Pct20','Pct30','Pct40','Pct50','Pct60','Pct70','Pct80','Pct90','Median','Location','southoutside')
        saveas(h1,sprintf('Results/Semi_Urb/%s Crit Branch Loads Case %d EV_penetration %d.png',time_of_year.name,selected_cases(c),EV_penetration));
        close(h1);
        hold off
        %save critical branch usefactor as csv
        writematrix(CriticalUsefactor(1:timesteps,iterations*(c-1)+1:iterations*c),sprintf('Results/Semi_Urb/%s Critical Branch Loads %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));


        h2=figure('visible','off');
        hold on
        plot(-mean(total_HH_load(:,iterations*(c-1)+1:iterations*c),2))
        plot(-mean(total_HH_load(:,iterations*(c-1)+1:iterations*c),2)-mean(total_EV_load(:,iterations*(c-1)+1:iterations*c),2))
        xlim([0 timesteps])
        xlabel('Time') 
        ylabel('Load (in MWh)')
        set(gca,'xtick',0:144:timesteps)
        set(gca,'xticklabel',{'Mon 0:00' 'Mon 12:00' 'Tue 0:00' 'Tue 12:00' 'Wed 0:00' 'Wed 12:00' 'Thu 0:00' 'Thu 12:00' 'Fri 0:00' 'Fri 12:00' 'Sat 0:00' 'Sat 12:00' 'Sun 0:00' 'Sun 12:00' })
        set(gca,'fontname','times', 'FontSize', 12)
        legend('Total Household Loads','Total Household and EV Loads','Location','southoutside')
        hold off
        saveas(h2,sprintf('Results/Semi_Urb/%s Loads Case %d EV Penetration %d.png',time_of_year.name, selected_case,EV_penetration));
        close(h2);
        if selected_case == 3 
            %save HH, EV as csv
            current_loads=[total_HH_load(:,iterations*(c-1)+1:iterations*c), total_HH_load(:,iterations*(c-1)+1:iterations*c) + total_EV_load(:,iterations*(c-1)+1:iterations*c)];
            writematrix(current_loads,sprintf('Results/Semi_Urb/%s HH Loads %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));
        end

        %plot transformer load
        h3=figure('visible','off');
        hold on
        data_plot = TransformerUsefactor(1:timesteps,iterations*(c-1)+1:iterations*c);
        plot(ones(timesteps))
        fanChart(1:size(data_plot,1), data_plot);
        ylim([0 1.5])
        xlim([0 timesteps])
        xlabel('Time') 
        ylabel('Load level (in %)')
        set(gca,'xtick',0:144:timesteps)
        set(gca,'xticklabel',{'Mon 0:00' 'Mon 12:00' 'Tue 0:00' 'Tue 12:00' 'Wed 0:00' 'Wed 12:00' 'Thu 0:00' 'Thu 12:00' 'Fri 0:00' 'Fri 12:00' 'Sat 0:00' 'Sat 12:00' 'Sun 0:00' 'Sun 12:00' })
        set(gca,'fontname','times', 'FontSize', 12)
        set(gcf,'position',[0,0,1200,1000])
        %legend('Pct10','Pct20','Pct30','Pct40','Pct50','Pct60','Pct70','Pct80','Pct90','Median','Location','southoutside')
        saveas(h3,sprintf('Results/Semi_Urb/%s Transformer Loads Case %d EV_penetration %d.png',time_of_year.name,selected_cases(c),EV_penetration));
        close(h3);
        %save critical transformer usefactor as csv
        writematrix(TransformerUsefactor(1:timesteps,iterations*(c-1)+1:iterations*c),sprintf('Results/Semi_Urb/%s Critical Transformer Loads %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));


        %plot min branch voltage
        h4=figure('visible','off');
        hold on
        data_plot = min_voltage(1:timesteps,iterations*(c-1)+1:iterations*c);
        plot(0.9*ones(timesteps))
        fanChart(1:size(data_plot,1), data_plot);
        ylim([0.8 1.05])
        xlim([0 timesteps])
        xlabel('Time') 
        ylabel('Voltage level (in %)')
        set(gca,'xtick',0:144:timesteps)
        set(gca,'xticklabel',{'Mon 0:00' 'Mon 12:00' 'Tue 0:00' 'Tue 12:00' 'Wed 0:00' 'Wed 12:00' 'Thu 0:00' 'Thu 12:00' 'Fri 0:00' 'Fri 12:00' 'Sat 0:00' 'Sat 12:00' 'Sun 0:00' 'Sun 12:00' })
        set(gca,'fontname','times', 'FontSize', 12)
        set(gcf,'position',[0,0,1200,1000])
        %legend('Pct10','Pct20','Pct30','Pct40','Pct50','Pct60','Pct70','Pct80','Pct90','Median','Location','southoutside')
        saveas(h4,sprintf('Results/Semi_Urb/%s Min voltage Case %d EV_penetration %d.png',time_of_year.name,selected_cases(c),EV_penetration));
        close(h4);
        close all;
        %save min voltage as csv
        writematrix(min_voltage(1:timesteps,iterations*(c-1)+1:iterations*c),sprintf('Results/Semi_Urb/%s Min voltage %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));


        %save critical branch IDs as csv
        writematrix(overloaded_branchID(:,iterations*(c-1)+1:iterations*c),sprintf('Results/Semi_Urb/%s Critical Branches %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));
    end
end
end



% POWER FLOW URBAN
% case selector:
% 1: charging loads induced by homecharging EVs
% 2: household electrical load profile
% 3: combined profiles of charging loads and households

            
urban.no_nodes = 108;

% specify loads per node (as specified in simBench, 0 for transformer buses)
loads_per_node = [3,0,0,0,1,0,2,2,0,3,1,3,3,2,2,2,3,1,1,1,1,2,4,2,3,1,1,2,4,1,4,3,1,3,3,1,3,1,2,3,2,1,1,1,2,2,2,2,2,3,3,2,1,1,2,2,2,2];
parfor time=1:number_of_seasons
    time_of_year = time_vars(time);
    [initial_EV_socs] = simulateSOCs(urban,time_of_year,housetype_charger_available, H16_car,car_segment_distribution, car_data, samples, speed,timestep_size, no_timesteps);
for EV_penetration = EV_pen
    min_voltage = zeros(timesteps,iterations*length(selected_cases));
    
    min_voltageID = zeros(timesteps,iterations*length(selected_cases));
    overloaded_branchID = zeros(timesteps,iterations*length(selected_cases));

    total_load = zeros(timesteps,iterations*length(selected_cases));
    total_HH_load = zeros(timesteps,iterations*length(selected_cases));
    total_EV_load = zeros(timesteps,iterations*length(selected_cases));
    branch_log = cell(timesteps,iterations*length(selected_cases));
    CriticalUsefactor=zeros(timesteps,iterations*length(selected_cases));
    CriticalBranchID=zeros(timesteps,iterations*length(selected_cases));
    TransformerUsefactor= zeros(timesteps,iterations*length(selected_cases));
    crit_fact_plot = nan(timesteps,iterations*length(selected_cases));
    hh_load_plot = nan(timesteps,iterations*11);
    hh_load_plot_idx = 1;
    for k = 1:iterations
        
        [nodes, charging_nodes] = assembleNodesN1(urban,time_of_year,EV_penetration,housetype_charger_available, H16_car,car_segment_distribution, car_data, samples, speed,timestep_size, no_timesteps, initial_EV_socs,house_scaling_factor);

        for l=1:length(selected_cases)
            selected_case=selected_cases(l);
            mpc = loadcase('SimBench_Urban6');
        %run load flow analysis without printing anything and enfore reactive limits
        mpopt = mpoption('pf.nr.max_it',100,'verbose',0,'out.all',0, 'pf.enforce_q_lims', 1);
        for i=1:timesteps
        % Pick #loads per node
        load_picker_idx = 1;
        % skip transformer busses of Urban6 Grid
        % skip bus 60 (n = 59)
        for n = 1:58
            % skip bus 9 (n = 9)
            if n ~= 9 %ok as no node attached
               % get profiles and convert in MW
               if loads_per_node(n) > 0
                   profile = nodes{n,selected_case}/1000;
               else
                   profile = zeros(timesteps,1);
               end
                total_HH_load(i,iterations*(l-1)+k) = total_HH_load(i,iterations*(l-1)+k) + nodes{n,2}(i)/1000;
                total_EV_load(i,iterations*(l-1)+k) = total_EV_load(i,iterations*(l-1)+k) + nodes{n,1}(i)/1000;
                 
                %Account for more then 1 load per node
                for more_nodes=1:loads_per_node(n)-1
                    profile = profile + nodes{load_picker_idx+more_nodes,selected_case}/1000;
                    total_HH_load(i,iterations*(l-1)+k) = total_HH_load(i,iterations*(l-1)+k) + nodes{load_picker_idx+more_nodes,2}(i)/1000;
                    total_EV_load(i,iterations*(l-1)+k) = total_EV_load(i,iterations*(l-1)+k) + nodes{load_picker_idx+more_nodes,1}(i)/1000;
                end
                load_picker_idx = load_picker_idx + loads_per_node(n);
                
               
                    % set bus
                    mpc.bus(n,2) = 1;
                    mpc.bus(n,3) = -profile(i)*0.93;
                    mpc.bus(n,4) = sqrt((profile(i))^2-(profile(i)*0.93)^2);
                    % set gen                
                    mpc.gen(n,2) = 0;
                    mpc.gen(n,9) = 0;                 
            end
        end
            result = runpf(mpc,mpopt);
            %Find critical branches
            Abranch=abs(complex(result.branch(2:end,14),result.branch(2:end,15))); %creates vector with apparent power values at the branches
            Capacity=mpc.branch(2:end,6); %creates a vector with the capity values of the branches
            Usefactor=Abranch./Capacity; %creates a vector with the usefactors at the branches
            CriticalUsefactor(i,iterations*(l-1)+k)=max(Usefactor); %determines the critical usefactor
            CriticalBranchID(i,iterations*(l-1)+k)=find(Usefactor == max(Usefactor(:)),1,'first'); %determines the critical branch
            
            if (result.branch(1,14)<0)
                TransformerUsefactor(i,iterations*(l-1)+k) = result.branch(1,14)./mpc.branch(1,6);
            else
                TransformerUsefactor(i,iterations*(l-1)+k) = abs(complex(result.branch(1,14),result.branch(1,15)))./mpc.branch(1,6);
            end
            %Min voltage violations were analyzed but did not occur,
            %hence not mentioned in final paper
            min_voltage(i,iterations*(l-1)+k) = min(result.bus(:, 8));
            
            min_voltageID(i,iterations*(l-1)+k) = find(result.bus(:, 8) == min(result.bus(:, 8)),1,'first');
            %only if usefactor >1 or min_voltage<0.9
            if CriticalUsefactor(i,iterations*(l-1)+k) > 1
                overloaded_branchID(i,iterations*(l-1)+k) = CriticalBranchID(i,iterations*(l-1)+k);
            elseif min_voltageID(i,iterations*(l-1)+k) < 0.9
                overloaded_branchID(i,iterations*(l-1)+k) = min_voltageID(i,iterations*(l-1)+k);
            elseif TransformerUsefactor(i,iterations*(l-1)+k) > 1
                overloaded_branchID(i,iterations*(l-1)+k) = 1;              
            end
           
            total_load(i,iterations*(l-1)+k) = sum(result.bus(:, 3));
            branch_log{i,1} = result.branch(:,:);  
        end
        end
        disp(k)
    end
     for c=1:length(selected_cases)
        selected_case=selected_cases(c);
        data_plot = CriticalUsefactor(1:timesteps,iterations*(c-1)+1:iterations*c);
        h1 = figure('visible','off');
        hold on
        plot(ones(timesteps))
        fanChart(1:size(data_plot,1), data_plot);
        ylim([0 1.5])
        xlim([0 timesteps])
        xlabel('Time') 
        ylabel('Load level (in %)')
        set(gca,'xtick',0:144:timesteps)
        set(gca,'xticklabel',{'Mon 0:00' 'Mon 12:00' 'Tue 0:00' 'Tue 12:00' 'Wed 0:00' 'Wed 12:00' 'Thu 0:00' 'Thu 12:00' 'Fri 0:00' 'Fri 12:00' 'Sat 0:00' 'Sat 12:00' 'Sun 0:00' 'Sun 12:00' })
        set(gca,'fontname','times', 'FontSize', 12)
        set(gcf,'position',[0,0,1200,1000])
        %legend('Pct10','Pct20','Pct30','Pct40','Pct50','Pct60','Pct70','Pct80','Pct90','Median','Location','southoutside')
        saveas(h1,sprintf('Results/Urban/%s Crit Branch Loads Case %d EV_penetration %d.png',time_of_year.name,selected_cases(c),EV_penetration));
        close(h1);
        hold off
        %save critical branch usefactor as csv
        writematrix(CriticalUsefactor(1:timesteps,iterations*(c-1)+1:iterations*c),sprintf('Results/Urban/%s Critical Branch Loads %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));


        h2=figure('visible','off');
        hold on
        plot(-mean(total_HH_load(:,iterations*(c-1)+1:iterations*c),2))
        plot(-mean(total_HH_load(:,iterations*(c-1)+1:iterations*c),2)-mean(total_EV_load(:,iterations*(c-1)+1:iterations*c),2))
        xlim([0 timesteps])
        xlabel('Time') 
        ylabel('Load (in MWh)')
        set(gca,'xtick',0:144:timesteps)
        set(gca,'xticklabel',{'Mon 0:00' 'Mon 12:00' 'Tue 0:00' 'Tue 12:00' 'Wed 0:00' 'Wed 12:00' 'Thu 0:00' 'Thu 12:00' 'Fri 0:00' 'Fri 12:00' 'Sat 0:00' 'Sat 12:00' 'Sun 0:00' 'Sun 12:00' })
        set(gca,'fontname','times', 'FontSize', 12)
        legend('Total Household Loads','Total Household and EV Loads','Location','southoutside')
        hold off
        saveas(h2,sprintf('Results/Urban/%s Loads Case %d EV Penetration %d.png',time_of_year.name, selected_case,EV_penetration));
        close(h2);
        if selected_case == 3 
            %save HH, EV as csv
            current_loads=[total_HH_load(:,iterations*(c-1)+1:iterations*c), total_HH_load(:,iterations*(c-1)+1:iterations*c) + total_EV_load(:,iterations*(c-1)+1:iterations*c)];
            writematrix(current_loads,sprintf('Results/Urban/%s HH Loads %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));
        end

        %plot transformer load
        h3=figure('visible','off');
        hold on
        data_plot = TransformerUsefactor(1:timesteps,iterations*(c-1)+1:iterations*c);
        plot(ones(timesteps))
        fanChart(1:size(data_plot,1), data_plot);
        ylim([0 1.5])
        xlim([0 timesteps])
        xlabel('Time') 
        ylabel('Load level (in %)')
        set(gca,'xtick',0:144:timesteps)
        set(gca,'xticklabel',{'Mon 0:00' 'Mon 12:00' 'Tue 0:00' 'Tue 12:00' 'Wed 0:00' 'Wed 12:00' 'Thu 0:00' 'Thu 12:00' 'Fri 0:00' 'Fri 12:00' 'Sat 0:00' 'Sat 12:00' 'Sun 0:00' 'Sun 12:00' })
        set(gca,'fontname','times', 'FontSize', 12)
        set(gcf,'position',[0,0,1200,1000])
        %legend('Pct10','Pct20','Pct30','Pct40','Pct50','Pct60','Pct70','Pct80','Pct90','Median','Location','southoutside')
        saveas(h3,sprintf('Results/Urban/%s Transformer Loads Case %d EV_penetration %d.png',time_of_year.name,selected_cases(c),EV_penetration));
        close(h3);
        %save critical transformer usefactor as csv
        writematrix(TransformerUsefactor(1:timesteps,iterations*(c-1)+1:iterations*c),sprintf('Results/Urban/%s Critical Transformer Loads %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));


        %plot min branch voltage
        h4=figure('visible','off');
        hold on
        data_plot = min_voltage(1:timesteps,iterations*(c-1)+1:iterations*c);
        plot(0.9*ones(timesteps))
        fanChart(1:size(data_plot,1), data_plot);
        ylim([0.8 1.05])
        xlim([0 timesteps])
        xlabel('Time') 
        ylabel('Voltage level (in %)')
        set(gca,'xtick',0:144:timesteps)
        set(gca,'xticklabel',{'Mon 0:00' 'Mon 12:00' 'Tue 0:00' 'Tue 12:00' 'Wed 0:00' 'Wed 12:00' 'Thu 0:00' 'Thu 12:00' 'Fri 0:00' 'Fri 12:00' 'Sat 0:00' 'Sat 12:00' 'Sun 0:00' 'Sun 12:00' })
        set(gca,'fontname','times', 'FontSize', 12)
        set(gcf,'position',[0,0,1200,1000])
        %legend('Pct10','Pct20','Pct30','Pct40','Pct50','Pct60','Pct70','Pct80','Pct90','Median','Location','southoutside')
        saveas(h4,sprintf('Results/Urban/%s Min voltage Case %d EV_penetration %d.png',time_of_year.name,selected_cases(c),EV_penetration));
        close(h4);
        close all;
        %save min voltage as csv
        writematrix(min_voltage(1:timesteps,iterations*(c-1)+1:iterations*c),sprintf('Results/Urban/%s Min voltage %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));


        %save critical branch IDs as csv
        writematrix(overloaded_branchID(:,iterations*(c-1)+1:iterations*c),sprintf('Results/Urban/%s Critical Branches %d EV_penetration %d.csv',time_of_year.name,selected_cases(c),EV_penetration));
    end
end
end




%%
%{ 
Bus Data Format
       1   bus number (positive integer)
       2   bus type
               PQ bus          = 1
               PV bus          = 2
               reference bus   = 3
               isolated bus    = 4
       3   Pd, real power demand (MW)
       4   Qd, reactive power demand (MVAr)
       5   Gs, shunt conductance (MW demanded at V = 1.0 p.u.)
       6   Bs, shunt susceptance (MVAr injected at V = 1.0 p.u.)
       7   area number, (positive integer)
       8   Vm, voltage magnitude (p.u.)
       9   Va, voltage angle (degrees)
   (-)     (bus name)
       10  baseKV, base voltage (kV)
       11  zone, loss zone (positive integer)
   (+) 12  maxVm, maximum voltage magnitude (p.u.)
   (+) 13  minVm, minimum voltage magnitude (p.u.)

 Generator Data Format
       1   bus number
   (-)     (machine identifier, 0-9, A-Z)
       2   Pg, real power output (MW)
       3   Qg, reactive power output (MVAr)       
       4   Qmax, maximum reactive power output (MVAr)
       5   Qmin, minimum reactive power output (MVAr)
       6   Vg, voltage magnitude setpoint (p.u.)
   (-)     (remote controlled bus index)
       7   mBase, total MVA base of this machine, defaults to baseMVA
   (-)     (machine impedance, p.u. on mBase)
   (-)     (step up transformer impedance, p.u. on mBase)
   (-)     (step up transformer off nominal turns ratio)
       8   status,  >  0 - machine in service
                    <= 0 - machine out of service
   (-)     (% of total VAr's to come from this gen in order to hold V at
               remote bus controlled by several generators)
       9   Pmax, maximum real power output (MW)
       10  Pmin, minimum real power output (MW)
   (2) 11  Pc1, lower real power output of PQ capability curve (MW)
   (2) 12  Pc2, upper real power output of PQ capability curve (MW)
   (2) 13  Qc1min, minimum reactive power output at Pc1 (MVAr)
   (2) 14  Qc1max, maximum reactive power output at Pc1 (MVAr)
   (2) 15  Qc2min, minimum reactive power output at Pc2 (MVAr)
   (2) 16  Qc2max, maximum reactive power output at Pc2 (MVAr)
   (2) 17  ramp rate for load following/AGC (MW/min)
   (2) 18  ramp rate for 10 minute reserves (MW)
   (2) 19  ramp rate for 30 minute reserves (MW)
   (2) 20  ramp rate for reactive power (2 sec timescale) (MVAr/min)
   (2) 21  APF, area participation factor

   Branch Data Format
       1   f, from bus number
       2   t, to bus number
   (-)     (circuit identifier)
       3   r, resistance (p.u.)
       4   x, reactance (p.u.)
       5   b, total line charging susceptance (p.u.)
       6   rateA, MVA rating A (long term rating)
       7   rateB, MVA rating B (short term rating)
       8   rateC, MVA rating C (emergency rating)
       9   ratio, transformer off nominal turns ratio ( = 0 for lines )
           (taps at 'from' bus, impedance at 'to' bus,
            i.e. if r = x = 0, then ratio = Vf / Vt)
       10  angle, transformer phase shift angle (degrees), positive => delay
   (-)     (Gf, shunt conductance at from bus p.u.)
   (-)     (Bf, shunt susceptance at from bus p.u.)
   (-)     (Gt, shunt conductance at to bus p.u.)
   (-)     (Bt, shunt susceptance at to bus p.u.)
       11  initial branch status, 1 - in service, 0 - out of service
   (2) 12  minimum angle difference, angle(Vf) - angle(Vt) (degrees)
   (2) 13  maximum angle difference, angle(Vf) - angle(Vt) (degrees)
           (The voltage angle difference is taken to be unbounded below
            if ANGMIN < -360 and unbounded above if ANGMAX > 360.
            If both parameters are zero, it is unconstrained.)
%}
