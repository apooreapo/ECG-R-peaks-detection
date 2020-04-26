function thresholdAdjustment

% Apostolou Orestis, for any questions contact me at
% orestisapostolou@yahoo.gr


% this function helps you to find the ideal threshold 
% it was used in order to reach to the choice of the 98th percentile
% must be in the same folder as fullEGCanalysis.m


    precision = zeros(5,15);
    recall = zeros(5,15);
    MSE = zeros(5,15);
    for i = 1:15
        j = 0.5+i./100;
        [precision(:,i), recall(:,i), MSE(:,i)] = ergasia2(j, 'rbio3.1');
    end
        
    [~,maxp] = max(precision');
    [~,maxr] = max(recall');
    [~,minm] = min(MSE');



end