function [initial_EV_socs,initial_BESS_socs] = simulateSOCs(agglomeration,time_of_year,pv_potential,BESS_data,housetype_charger_available, H16_car,car_segment_distribution, car_data, samples, speed,timestep_size, no_timesteps)
%SIMULATE SOCS for initialization
%  
% set PV pen to 100% to ensure that every house has a BESS
pv_pen.single_house = 1;
pv_pen.double_house= 1;
pv_pen.multi_house= 1;
load_nodes = assembleNodes3(agglomeration,time_of_year,pv_potential,pv_pen,BESS_data,100,housetype_charger_available, H16_car,car_segment_distribution, car_data, samples, speed,timestep_size, no_timesteps,0.5,0.5,0);
for i=1:length(load_nodes)
    car = EV(timestep_size, 60, 50, 0.2, speed,time_of_year.consumption_factor, 11, 0.9);
    simulateEV(car,samples,11,1,5,no_timesteps);
    initial_EV_socs(i) = car.returnSOC();
    storage = BESS(BESS_data.capacity,BESS_data.SOC,BESS_data.charging_efficiency,timestep_size);
    storage.simulate(no_timesteps,load_nodes{i,5});
    initial_BESS_socs(i) = storage.returnSOC();
end
end

