classdef BESS < handle
    %CAR Summary of this class goes here
    %
    properties
        batterycap
        soc
        charging_eff
        timefactor
        timestep_size
    end
    
    methods
        function obj = BESS(batterycap, soc, charging_eff, timestep_size)
            % Construct an instance of this class
            %    :param batterycap: the usable batterycapacity of the car in kWh
            %    :param soc: the initial state of charge of the battery in % (0-100)
            %    :param charging_eff: factor to account for energy losses during charging, default = 1   
            %    :param timestep_size: the size of each timestep in min
            obj.timestep_size = timestep_size;
            obj.timefactor = (timestep_size/60); %because values are in kWh
            obj.batterycap = batterycap;
            obj.soc = soc;
            obj.charging_eff = charging_eff;
        end 

        function [power, power_interaction] = simulate(obj, no_timesteps, load_profile)
            if length(load_profile)~=no_timesteps
                warning("Load profile does not match the number of timesteps");
            end
            temp_power = zeros(no_timesteps,1);
            for i = 1:no_timesteps
                if load_profile(i) >= 0
                    temp_power(i) = obj.charge(load_profile(i));
                else
                    temp_power(i) = obj.discharge((-1)*load_profile(i));
                end
            end
            
            % return drawn/provided power; negative = sink, positive = supply
            power_interaction = -temp_power;
            power = load_profile - temp_power;
        end

        function [power] = discharge(obj, discharge_power)
            %discharge Summary of this method goes here
            % simulates driving for one timestep
            % returns soc, power consumption
            consumed_energy = obj.timefactor * discharge_power;
            obj.soc = obj.soc - 100*(consumed_energy/obj.batterycap);
            if obj.soc<=0
                obj.soc = 0;
                power = 0;
            else
                power = -discharge_power;
            end
        end
        
        function [power] = charge(obj, charging_power)
            %charge Summary of this method goes here
            %simulates charging for one timestep
            %returns soc, power_consumption                
            charged_power = charging_power * obj.charging_eff * obj.timefactor;
            obj.soc = obj.soc + 100*(charged_power/obj.batterycap);
            power = charging_power;
            if obj.soc>=100
                obj.soc = 100;
                power = 0;
            end          
        end

        function soc = returnSOC(obj)
            soc = obj.soc;
        end
    end
end

