function [newVals] = boxCarAvg(data, reductionFactor) % Give this function your data as a col vector please :)
    newVals = zeros(length(data)/reductionFactor,1); %preallocate size of new values
    numIters = 1; % iteration count
    for i = 1:reductionFactor:length(data)-reductionFactor % takes step of size of reduction factor 
            selectedVals = data(i:i+reductionFactor-1,1); %select reductionFactor number of values
            newVals(numIters,1) = sum(selectedVals)/length(selectedVals); % take the average of selected values and assign to newVals
            numIters = numIters+1; % increment iterations
    end
end