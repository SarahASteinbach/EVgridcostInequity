classdef EV < handle
    %CAR Summary of this class goes here
    %
    properties
        timestep_size
        timefactor
        batterycap
        soc
        energy_consumption
        speed
        consumption_factor
        charging_power
        charging_eff
    end
    
    methods
        function obj = EV(timestep_size, batterycap, soc, energy_consumption, speed,consumption_factor, charging_power, charging_eff)
            %CAR Construct an instance of this class
            %    :param name: name of the car model
            %    :param timestep_size: the size of each timestep in min
            %    :param batterycap: the usable batterycapacity of the car in kWh
            %    :param soc: the initial state of charge of the battery in % (0-100)
            %    :param energy_consumption: the average energy_consumption of the car in kWh/km
            %    :param avg_speed: the average speed of the car in km/h
            %    :param consumption_factor: consumption factor depending on ambient_temp, environment, driving style
            %    :param charging_power: charging power of the car in kW
            %    :param charging_eff: factor to account for energy losses during charging 
            obj.timestep_size = timestep_size;
            obj.timefactor = (timestep_size/60);
            obj.batterycap = batterycap;
            obj.soc = soc;
            obj.energy_consumption = energy_consumption;
            obj.speed = speed;
            obj.consumption_factor = consumption_factor;
            obj.charging_power = charging_power;
            obj.charging_eff = charging_eff;
        end     
        

        function [soc, power] = drive(obj,timepoint)
            %drive Summary of this method goes here
            % simulates driving for one timestep
            % returns soc, power consumption
            drawn_speed = obj.drawSpeed(timepoint);
            consumed_energy = drawn_speed * obj.timefactor * obj.energy_consumption * obj.consumption_factor;
            obj.soc = obj.soc - 100*(consumed_energy/obj.batterycap);
            obj.selfDischarge();
            soc = obj.soc;
            power = 0;
        end
        
        function [soc, power] = charge(obj)
            %charge Summary of this method goes here
            %simulates charging for one timestep
            %returns soc, power_consumption                
            charged_power = obj.charging_power * obj.charging_eff * obj.timefactor;
            obj.soc = obj.soc + 100*(charged_power/obj.batterycap);
            obj.selfDischarge();
            power = obj.charging_power;
            if obj.soc>=100
                obj.soc = 100;
                power = 0;
            end
            soc = obj.soc;            
        end

        function [soc, power] = park(obj)
            %park Summary of this method goes here
            %simulates parking for one timestep
            %returns soc, power_consumption 
            obj.selfDischarge();
            soc = obj.soc;
            power = 0;
        end

        function none = selfDischarge(obj)
            obj.soc = obj.soc - (obj.timefactor/24);
        end
        

        function speed = drawSpeed(obj,timepoint)
            day = ceil (timepoint/(1440/obj.timestep_size));
            daytime = timepoint - (day-1)*(1440/obj.timestep_size);
            if day <= 5 
                if ((360/obj.timestep_size) < daytime) && (daytime <= (720/obj.timestep_size))
                    speed= randsample(obj.speed.morning_wd.SPEED,1);
                elseif ((720/obj.timestep_size) < daytime) && (daytime <= (1200/obj.timestep_size))
                    speed= randsample(obj.speed.afternoon_wd.SPEED,1);
                else
                    speed= randsample(obj.speed.night_wd.SPEED,1);
                end
            else
                if ((360/obj.timestep_size) < daytime) && (daytime <= (720/obj.timestep_size))
                    speed= randsample(obj.speed.morning_we.SPEED,1);
                elseif ((720/obj.timestep_size) < daytime) && (daytime <= (1200/obj.timestep_size))
                    speed= randsample(obj.speed.afternoon_we.SPEED,1);
                else
                    speed= randsample(obj.speed.night_we.SPEED,1);
                end
            end
        end

        function soc = returnSOC(obj)
            soc = obj.soc;
        end

        function none = setConsumptionFactor(obj, new_consumption_factor)
            obj.consumption_factor = new_consumption_factor;
        end

        function none = setChargingPower(obj, new_charging_power)
            obj.charging_power = new_charging_power;
        end
    end
end

