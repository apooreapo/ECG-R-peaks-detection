function [TP, FP, FN, precision, recall, matchMat] = evaluation(f1,n1,n2,margin)


% f1 is a matrix with 2 columns. The first one has the locations of the
% real annotations, and the second one has the predicted annotations. 
% n1 is the size of the real annotations and n2 is the size of the
% predicted ones
% if n1 is not equal to n2, the column with the least elements must be
% zero padded in the ending. Both columns must be sorted increasingly.
% Margin is the max margin for which a distance between predicted and
% actual annotation is considered as true prediction

% TP : True Positive
% FP : False Positive
% FN : False Negative
% matchMat : The matrix with the matching of real and predicted annotations
    
    [N,~] = size(f1);
    matchMat = zeros(N,2);
    count = 1;
    j = 1;
    i=1;
    while (i<= n1 && j <=n2)
        if (abs(f1(i,1) - f1(j,2)) <= margin)
            matchMat(count,1) = i;
            matchMat(count,2) = j;
            count = count+1;
            i = i+1;
            j = j+1;
        else
            if f1(i,1)<f1(j,2)
                i = i+1;
            else
                j = j+1;
            end
        end
    end
    matchMat( ~any(matchMat,2), : ) = [];   % delete zero rows
    [TP,~] = size(matchMat);
    FN = n1 - TP;
    FP = n2 - TP;
    precision = TP / (TP + FP);
    recall = TP / (TP + FN);
    
end

function MSE = calculateMSE(signal, matchMat)

% signal is the ecg we are interested in

    [n1,~] = size(matchMat);
    MSE = 0;
    for i = 1: n1
        MSE = MSE + (signal(matchMat(i,1)) - signal(matchMat(i,2))).^2;
    end
    MSE = MSE./n1;


end