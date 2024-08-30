function [chain,idx,total_trip_duration] = generateDailyTrips(chain,wotag,first_trips,transitionMatrixNew,parkingTimes,no_trips,dauer_prob,prev_end_idx)
%DAILYTRIPS Summary of this function goes here
%   Detailed explanation goes here
%disp("Generate Trips WOTAG " + wotag + " wotag idx: " + weekdayIndex(wotag) + "     | no_trips: " +  no_trips + " | prev end idx: " + prev_end_idx)
% make sure to start the first trip after the last trip of the prev day ended
%Creates chain from first trip of the day until the last trip, idx
%represents arrival time after a trip was made (in output of last trip of
%the day)
while 1
    starttime = randsample(first_trips.abzeit_i_ts,1)+weekdayIndex(wotag);
    if starttime > prev_end_idx
        break
    end
end
% select the parking times based on daily trip frequencies and aggregate frequencies above 7 to maintain stable results
if no_trips > 7     
    parkingTimes = parkingTimes(parkingTimes(:,4)>7,:);
else 
    parkingTimes = parkingTimes(parkingTimes(:,4)==no_trips,:);
end
% assign home for the whole day if no trips conducted on that day
total_trip_duration = 0;
idx = prev_end_idx;
%disp("No trips" + no_trips)
if no_trips == 0
    chain(prev_end_idx+1:weekdayIndex(wotag+1)) = 1;
else if no_trips == 1   
    %empirically sample trip duration from data
    tripDuration = randsample(first_trips.duration_ts,1);
    chain(prev_end_idx+1:starttime) = 1;
    chain(starttime+1:starttime+tripDuration) = 0;
    idx = starttime+tripDuration;
    total_trip_duration = total_trip_duration + tripDuration;
else
    %iterate through the trips
    for tripidx = 1:no_trips
        if idx > 2160
            break;
        end
        %first trip of day starts at home
        if tripidx == 1
            tripDuration = randsample(first_trips.duration_ts,1);
            nextState = drawNextState(starttime,1,transitionMatrixNew);
            %assign home state until the starttime and the driving state for the trip time
            chain(prev_end_idx+1:starttime) = 1;
            chain(starttime+1:starttime+tripDuration) = 0; 
            %disp("Trip duration should: " +  tripDuration + " | trip duration is: " + (starttime+tripDuration - starttime+1-1))
            total_trip_duration = total_trip_duration + tripDuration;
            % if only two trips on that day make sure that its a round trip and that duration equals
            if no_trips == 2
                parkingDuration = randsample(parkingTimes(parkingTimes(:,1)==nextState,2),1);
                chain(starttime+tripDuration+1:starttime+tripDuration+parkingDuration) = nextState;             
                chain(starttime+tripDuration+parkingDuration+1:starttime+tripDuration+parkingDuration+tripDuration) = 0;
                idx = starttime+tripDuration+parkingDuration+tripDuration;
                %disp("Trip duration should: " +  tripDuration + " | trip duration is: " + (starttime+tripDuration+parkingDuration+1-starttime+tripDuration+parkingDuration+tripDuration-1))
                total_trip_duration = total_trip_duration + tripDuration;
                break
            % if more then 2 trips: draw parkingDuration based on the state and assign the nextState for the whole parkingDuration
            else
                parkingDuration = randsample(parkingTimes(parkingTimes(:,1)==nextState,2),1);
                chain(starttime+tripDuration+1:starttime+tripDuration+parkingDuration)= nextState;
                idx = starttime+tripDuration+parkingDuration;
            end
        %if not the last trip
        elseif tripidx ~= no_trips
            %determine next state based on current state and time
            % mod avoids index out of bounds for transitionMatrix
            nextState = drawNextState(mod(idx-1,2015)+1,chain(idx),transitionMatrixNew);
            tripDuration = myrandsrc(1,1,[1:26; cumsum(dauer_prob)]); 
            parkingDuration = randsample(parkingTimes(parkingTimes(:,1)==nextState,2),1);
            chain(idx+1:idx+tripDuration) = 0;
            chain(idx+1+tripDuration:idx+tripDuration+parkingDuration) = nextState;
            idx = idx+tripDuration+parkingDuration;
            %disp("Trip duration should: " +  tripDuration + " | trip duration is: " + (idx+1-idx+tripDuration-1))
            total_trip_duration = total_trip_duration + tripDuration;
        %last trip of the day ends at home
        else
            %adjusted
            tripDuration = myrandsrc(1,1,[1:26; cumsum(dauer_prob)]); 
            chain(idx+1:idx+tripDuration) = 0;
            idx = idx+tripDuration;
            %disp("Trip duration should: " +  tripDuration + " | trip duration is: " + (idx+1-idx+tripDuration-1))
            total_trip_duration = total_trip_duration + tripDuration;
        end  
        %disp("Trip_idx: " + tripidx + " | current duration: " + total_trip_duration)
    end
end
end
