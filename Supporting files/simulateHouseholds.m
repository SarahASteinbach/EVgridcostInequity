function [simulated_households] = simulateHouseholds(no_simulated_households,agglomeration, time_of_year, H16, car_segment_distribution, car_data, speed,timestep_size,no_timesteps,initial_EV_socs)
%Simulates electrical loads for a household for a duration of 1 week
% RETURN
% row: individual household
% col 1: number of people in household
% col 2: number of cars in that household
% col 3: assigned vehicles in that household
% col 4: household electrical load profile
simulated_households = cell(no_simulated_households,4);

% Sample a householdsize
household_size = myrandsrc(no_simulated_households,1,[1 2 3 4 5; cumsum(agglomeration.pers_per_hh)]);
for i = 1:no_simulated_households
    
    % Sample if car in household
    if_car=binornd(1,agglomeration.hh_with_car);
    % Sample no of cars of this household
    no_cars = if_car*randsample(H16(H16.HHGRO==household_size(i),:).PKWHH,1);

    % Sample load profile for household depending on size and time of the
    % year (from load profile Generator)
    household_profiles = time_of_year.household_profiles(time_of_year.household_profiles(:,1)==household_size(i),2:end).';
    household_load = household_profiles(:,randsample(1:width(household_profiles),1));
    % create 5 minute timestep data
    household_load = sum(reshape(household_load,timestep_size,no_timesteps))*(60/timestep_size);
    % sample car models
    %this section is not executed if no_cars=0
    carTypes = myrandsrc(no_cars,1,[2 3 4 5 6; cumsum(car_segment_distribution)]);
    cars = cell(1,no_cars);
    for j = 1:no_cars
        batterycap = car_data{1,carTypes(j)};
        energy_consumption = car_data{2,carTypes(j)};
        % random draw SOC
        soc = randsample(initial_EV_socs,1);
        cars{1,j} = EV(timestep_size, batterycap, soc, energy_consumption, speed,time_of_year.consumption_factor, 11, 0.9);
    end
    simulated_households{i,1} = household_size(i);
    simulated_households{i,2} = no_cars;
    simulated_households{i,3} = cars;
    simulated_households{i,4} = -household_load;
end
end
