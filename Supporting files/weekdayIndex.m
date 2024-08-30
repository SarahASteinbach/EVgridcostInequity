function [idx] = weekdayIndex(WOTAG)
%returns the index of the weekday
    idx = 288*(WOTAG-1)+1;
end

