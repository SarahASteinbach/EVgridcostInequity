function [nextState] = drawNextState(timeidx,currentState,transitionMatrixNew)
%DRAWNEXTSTATE Summary of this function goes here
%   Detailed explanation goes here
    %draw transitionMatrix next state based on state at the current timestep
    x = transitionMatrixNew(timeidx,currentState,:);
    x = x(:).';
    %draw next State based on transitonMatrix
    nextState = myrandsrc(1,1,[1 2 3; cumsum(x)]);
end


