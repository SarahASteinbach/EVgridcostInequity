function [status] = returnStatus(zweck,prev_status)
%% encodes the status based on trip purpose
% 0=Driving (D)
% 1=Home (H)
% 2=Work (W)
% 3=Other (O)
% 4=OtherTransport (not needed in current version)
    % assign default for prev_status
    if nargin == 1
        prev_status = 1;
    end
    if zweck == 1 || zweck == 3
        status = 2;
    elseif zweck == 7 || zweck == 10
        status = 1;
    elseif zweck == 77
        status = prev_status;
    else
        status = 3;
    end
end