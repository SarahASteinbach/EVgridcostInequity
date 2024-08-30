function [states] = assignStatus(startindex, endindex, weekday_index, UUID, status, states)
%% assigns the status to the state variable between startIdx and endIdx
%handle case where new week beginns and assign status to all indexes between start and end
    %disp("UUID: " + UUID + " | startindex: " + startindex + " | endindex: " + endindex)
    if endindex < startindex
        states{startindex:2016, UUID} = status;
        states{1:endindex, UUID} = status;
    else
        states{startindex:endindex, UUID} = status;
    end
end