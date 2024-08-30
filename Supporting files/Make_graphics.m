%Example graph for Rural 40/40 load visual

Plot_Loads = table2array(readtable('Y:\TUM-PC\Modelling Sarah\Results\For Paper\Dif Penetration\Poor\Rural'));
Plot_HH = -mean(Plot_Loads(:,1:iterations),2);
Plot_HH_EV = -mean(Plot_Loads(:, iterations+1:iterations*2),2);
Plot_PV = mean(Plot_Loads(:,iterations*2+1:iterations*3),2);
Plot_BESS = mean(Plot_Loads(:,iterations*3+1:iterations*4),2);

timesteps=288;

h1=figure('visible','off');
hold on
plot(Plot_HH,'Color',' #757575','LineStyle','--')
plot(Plot_HH_EV,'Color',' #757575 ')
plot(Plot_HH_EV-Plot_PV)
plot(Plot_HH_EV-Plot_PV-Plot_BESS,'Color',' blue ')
plot(zeros(timesteps))
xlim([0 timesteps])
xlabel('Time') 
ylabel('Load (in MWh)')
legend('Total HH Loads','Total HH and EV Loads','HH, EV loads and PV generation','HH, EV loads, PV generation and BESS in-/outflow', 'Location','southoutside')
hold off
saveas(h1,sprintf('Results/For paper/Graphics/Loadexample.png'));
close(h1);



timesteps=2016;

h2=figure('visible','off');
hold on
plot(Plot_HH_EV,'Color',' #757575 ')
xlim([0 timesteps])
xlabel('Time') 
ylabel('Load (in MWh)')
legend('Total HH and EV Loads', 'Location','southoutside')
hold off
saveas(h2,sprintf('Results/For paper/Graphics/Loadexample_week.png'));
close(h2);