function [prec, rec, MSEfull] = fullECGanalysis(thres,wavelet)

% Apostolou Orestis, for any questions contact me at
% orestisapostolou@yahoo.gr

% R PEAK DETECTION ALGORITHM
%
% ergasia2.m must be in the same folder as 106.dat, 107.dat ...

% read the ecgs and their annotations

% run as ergasia2(), or ergasia2(thres), where thres must be between 0 and
% 1, or as ergasia2(thres, wavelet) where wavelet is the wavelet you want
% to use. If you use thres = 0 the algorithm finds automatically the
% threshold from the 98th percentile

% By default, the output of this script is the precision of each ECG R peak
% prediction. If you want to see other evaluation metrics, you either need
% to put a breakpoint at the end of the script, or modify it appropriately

if nargin == 0
    thres = 0;
    wavelet = 'fk4';
elseif nargin == 1
    wavelet = 'fk4';
end

ann106=rdann('106.dat','atr',1);
[ecg106, Fs, tm] = rdsamp('106.dat',1);

ann107=rdann('107.dat','atr',1);
[ecg107, Fs, tm] = rdsamp('107.dat',1);

ann113=rdann('113.dat','atr',1);
[ecg113, Fs, tm] = rdsamp('113.dat',1);

ann202=rdann('202.dat','atr',1);
[ecg202, Fs, tm] = rdsamp('202.dat',1);

ann232=rdann('232.dat','atr',1);
[ecg232, Fs, tm] = rdsamp('232.dat',1);

prec = zeros(5,1);
rec = zeros(5,1);
MSEfull = zeros(5,1);

[prec(1), rec(1), MSEfull(1)] = analyzeECG(ecg106, ann106, thres, wavelet);
[prec(2), rec(2), MSEfull(2)] = analyzeECG(ecg107, ann107, thres, wavelet);
[prec(3), rec(3), MSEfull(3)] = analyzeECG(ecg113, ann113, thres, wavelet);
[prec(4), rec(4), MSEfull(4)] = analyzeECG(ecg202, ann202, thres, wavelet);
[prec(5), rec(5), MSEfull(5)] = analyzeECG(ecg232, ann232, thres, wavelet);

%figure
%plot(prec);

%figure
%plot(rec);

%figure
%plot(MSEfull);

mPrecision = mean(prec);
mRecall = mean(rec);
mMSE = mean(MSEfull);


% Place a breakpoint here to view precision, recall and F1 measures fot the
% 5 ecgs in the variables menu


end



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

function MSE = calculateMSE(f1, matchMat)

    [n1,~] = size(matchMat);
    MSE = 0;
    for i = 1: n1
        MSE = MSE + (f1(matchMat(i,1),1) - f1(matchMat(i,2),2)).^2;
    end
    MSE = MSE./n1;
    MSE = sqrt(MSE);


end

function [precision, recall, MSE] = analyzeECG(ecg, ann, thres, wavelet)
    
    N = length(ecg);
    x = 1:1:N;
    y(1:N) = thres; % threshold
    %figure
    %plot(x,ecg);
    %hold on
    %plot(x(ann(:)),ecg(ann(:)),'r*');

    %[cA1,cD1] = dwt(ecg106,'coif1');
    %A1coif = upcoef('a',cA1,'coif1',1,N);
    %D1 = upcoef('d',cD1,'coif1',1,N);

    % wavelet decomposition into 2 levels, 'fk4' is the wavelet family

    [C,L] = wavedec(ecg,2,wavelet);
    [cD1, cD2] = detcoef(C,L,[1,2]);

    % D1 is first detail (higher freqs) and D2 is the second

    D1 = wrcoef('d',C,L, wavelet,1);
    D2 = wrcoef('d',C,L, wavelet,2);

    %figure
    %plot(ecg(1:10000));
    %hold on
    %plot(D2(1:10000));

    % we use D2 because it has lower noise

    mx = max(D2);
    mn = min(D2);
    D2normal(:) = D2(:)-mn;
    D2normal(:) = D2normal(:)/(mx-mn);
    if thres == 0
        thres = prctile(D2normal,98);
    end

    [peaks, locs] = findpeaks(D2normal,'MinPeakDistance',90,'MinPeakHeight',thres);

    % 90 samples as minPeakDistance means that this application ignores pulses
    % with more than 240 bpm

    n1 = length(locs);
    n2 = length(ann);
    f1 = zeros(max(n1,n2),2);
    locsT = locs';

    locsFinal = zeros(size(locsT));
    for i = 1:n1
         limL = max(locsT(i)-7,1);
         limU = min(locsT(i)+7,N);
         [~,locsFinal(i)] = max(ecg(limL:limU));
         locsFinal(i) = locsFinal(i)+limL-1;
    end




    if n1>n2
        ann(end+1:n1,1)=0;
    else
        locsFinal(end+1:n2,1) = 0;
    end


    f1(:,1) = ann;
    f1(:,2) = locsFinal;

    % f1 first column is the real R annotations, and second column is the
    % predicted ones

    locsFinal = nonzeros(locsFinal);
    %figure
    %title("o is real, x is predicted");
    %hold on
    %plot(x,ecg);
    %hold on
    %plot(x,D2);
    %hold on
    %plot(x,y);
    %hold on
    %plot(x(locsFinal(:)),ecg(locsFinal(:)),'r*');
    %hold on
    %plot(x(ann(:)),ecg(ann(:)),'ko', 'Color', [0 0 0]);

    [TP, FP, FN, precision, recall, matchMat] = evaluation(f1,n1,n2,4);
    MSE = calculateMSE(f1, matchMat); 



end