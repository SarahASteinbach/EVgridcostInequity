function [powerconsumption, socs] = simulateEV(EV, samples, charging_power, decision_factor, timestep_size, no_timesteps)
%create the energy consumption profiles of the EVs of a node

    powerconsumption=nan(no_timesteps,1);
    socs=zeros(no_timesteps,1);
    % draw a random sample
    k = floor((width(samples)-1).*rand(1,1) + 1);
    EV.setChargingPower(charging_power)
    % a sets the coefficient for the charging decision
    a = 0.15;
    factor = decision_factor;
    while (any(socs(:) <= 0))
        chargingDecision = false;
        decisionMade = true;
        for i=1:no_timesteps
            if samples(i,k)==0
                [soc, power] = EV.drive(i);
                chargingDecision = false;
                decisionMade = false;
            elseif samples(i,k)==1 
                %decide if charging should be done
                if ~decisionMade
                    %determine duration of the parking stay
                    for endParkingidx=i:no_timesteps
                        if samples(endParkingidx,k)~=1
                            break
                        end
                    end
                    % only consider charging if parking duration > 10 min
                    if endParkingidx-i > 2
                        current_soc = EV.returnSOC();                        
                        pcharge = min(1,1-factor/(1+exp(-a*(current_soc-60))));
                        % decide if charging
                        if rand < pcharge
                            chargingDecision = true;
                        end       
                    end
                    decisionMade = true;
                end
                if chargingDecision
                    [soc, power] = EV.charge();
                else
                    [soc, power] = EV.park();
                end
             
            elseif samples(i,k)==2
                [soc, power] = EV.park();
                chargingDecision = false;
                decisionMade = false;
       
            else 
                [soc, power] = EV.park();
                chargingDecision = false;
                decisionMade = false;
            end
            powerconsumption(i) = power;
            socs(i) = soc;
        end  
    end
end 


