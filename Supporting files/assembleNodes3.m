function [nodes,charging_nodes, house_loads] = assembleNodes3(agglomeration,time_of_year,pv_potential,pv_penetration, BESS_data,EV_penetration,charger_available,H16,car_segment_distribution, car_data, samples, speed, timestep_size, no_timesteps,initial_EV_socs,initial_BESS_socs,house_scaling_factor)
%ASSEMBLENODES 
% Params
% agglomeration: the type of environment (rural, ...)
% time_of_year: e.g. march
% pv_potential: the housesize dependent pv potential
% BESS_data: Data to init the BESS
% EV_pentration: The EV Penetration in % (e.g. 30)
% charger_available: Data specifying the availability of charging stations
% H16_car: Dataset containing household size and no of cars
% car_segment_distribution: Specifies the segments of EVs
% car_data: Data to init the vehicles
% samples: State Chain to determine the vehicles location
% speed: Dataset containing speed samples
% RETURN
% of nodes:
% col 1: charging loads induced by charging EVs
% col 2: household electrical load profile
% col 3: combined profiles of charging loads and households
% col 4: pv generation profile for the node
% col 5: combined electrical profiles (positive values = generation, negative values = load)
% col 6: BESS interaction drawn/provided power
% col 7: electrical profile with BESS interaction
%and charging nodes
no_nodes = agglomeration.no_nodes;
nodes = cell(no_nodes,7);
charging_nodes = cell(no_nodes,1);
% itermediary storage for loads on house level
house_loads = cell(no_nodes,2,7);
%Add an additional house to the node in x% of cases to match load
%profiles
if house_scaling_factor>=0 
    no_houses = 1 + binornd(1,house_scaling_factor,1,no_nodes);
else
    no_houses = 1 - binornd(1,-house_scaling_factor,1,no_nodes);
end
%disp("No houses: " + sum(no_houses))
% sample number of households per node and simulate required amount of households
no_households = myrandsrc(sum(no_houses),1,[agglomeration.hh_per_house_amount; cumsum(agglomeration.hh_per_house)]);
simulated_households = simulateHouseholds(sum(no_households),agglomeration,time_of_year,H16,car_segment_distribution, car_data, speed, timestep_size, no_timesteps,initial_EV_socs);

households_counter = 1;
house_counter=1;

no_pvs=0;
for i = 1:no_nodes
    % scale PV production according to available pv capacity for the
    % house and set likelihood of home charger availability
    for k=1:no_houses(i)
        house_with_ev = false;
        pv_factor=0;
        if no_households(house_counter) == 1
            home_charging = agglomeration.private_parking_prob * charger_available.single_house;
        elseif no_households(house_counter) == 2
            home_charging = agglomeration.private_parking_prob * charger_available.double_house;
        else      
            home_charging = agglomeration.private_parking_prob * charger_available.multi_house;
        end
        % decide if charging can be done at this house
        if rand <= home_charging
            home_charging = true;
            charging_power = 11;
            decision_factor = 1;
        else
            home_charging = false;
            charging_power = 22;
            decision_factor = 1.2;
        end
        %household level
        for j = 1:no_households(house_counter)
            %if the household owns EVs simulate EV
            no_cars = simulated_households{households_counter,2};
            cars = simulated_households{households_counter,3};
            no_EVs = binornd(no_cars,EV_penetration/100);
            powerconsumption = zeros(no_timesteps,1);
            for l=1:no_EVs
                    house_with_ev = true;
                    [simulatedPower, returnedSoc] = simulateEV(cars{l},samples,charging_power,decision_factor,timestep_size,no_timesteps);
                    powerconsumption = powerconsumption + simulatedPower;
            end
             %Let charge while at home
            home_charging = 1;
            if home_charging
                if isempty(house_loads{i,k,1})
                    house_loads{i,k,1} = -powerconsumption;
                else
                    house_loads{i,k,1} = house_loads{i,k,1} - powerconsumption;
                end
                charging_nodes{i,1} = zeros(no_timesteps,1);
            else
                house_loads{i,k,1} = zeros(no_timesteps,1);              
                if isempty(charging_nodes{i,1})
                    charging_nodes{i,1} = -powerconsumption;
                else
                    charging_nodes{i,1} = charging_nodes{i,1} - powerconsumption;
                end
            end     
            if isempty(house_loads{i,k,2})
                house_loads{i,k,2} = simulated_households{households_counter,4}.';
            else
                house_loads{i,k,2} = house_loads{i,k,2} + simulated_households{households_counter,4}.';
            end
            households_counter = households_counter + 1;
        end
        %end household level

        % PV on house level
        % 31% more likely to have PV in case of EV
        if house_with_ev
            %PV_penetration
            pv_house_factor = 1.31;
        else
            pv_house_factor = (1-EV_penetration/100*1.31)/(1-EV_penetration/100);
        end
        if no_households(house_counter) == 1
            if rand < pv_penetration.single_house*pv_house_factor
                pv_factor = pv_potential.single_house/10;
                no_pvs=no_pvs+1;
            end
        elseif no_households(house_counter) == 2
            if rand < pv_penetration.double_house*pv_house_factor
                pv_factor = pv_potential.double_house/10;
                no_pvs=no_pvs+1;
            end
        else
            if rand < pv_penetration.multi_house*pv_house_factor
                pv_factor = pv_potential.multi_house/10;
                no_pvs=no_pvs+1;
            end        
        end
        %if no EVs set demand to 0
        if isempty(house_loads{i,k,1})
            house_loads{i,k,1} = zeros(no_timesteps,1);
        end
        house_loads{i,k,3} = house_loads{i,k,1} + house_loads{i,k,2};
        %assign PV gen capacity for the house
        house_loads{i,k,4} = time_of_year.pv.' * pv_factor;
        house_loads{i,k,5} = house_loads{i,k,3} + house_loads{i,k,4};

        % only battery capacity if solar PV exists
        if pv_factor > 0
            storage = BESS(BESS_data.capacity, randsample(initial_BESS_socs,1),BESS_data.charging_efficiency,timestep_size);
            [house_loads{i,k,7}, house_loads{i,k,6}] = storage.simulate(no_timesteps,house_loads{i,k,5});
        else
            house_loads{i,k,6} = zeros(no_timesteps,1);
            house_loads{i,k,7} = house_loads{i,k,5};
        end
        house_counter = house_counter + 1;
    end
    if isempty(charging_nodes{i,1})
        charging_nodes{i,1} = zeros(no_timesteps,1);
    end
end
% create node profiles by summing up over all houses at a node
nodes = cell(no_nodes,7);
for i = 1:no_nodes
    if no_houses(i)==0
        for j = 1:7
            nodes{i,j} = zeros(no_timesteps,1);
        end
    else
        for k=1:no_houses(i)
            for j = 1:7
                if isempty(nodes{i,j})
                    nodes{i,j} = house_loads{i,k,j};
                else
                    nodes{i,j} = nodes{i,j} + house_loads{i,k,j};
                end
            end
        end
    end
end

end


