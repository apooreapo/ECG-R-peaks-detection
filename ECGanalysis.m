function ECGanalysis

% Apostolou Orestis, for any questions contact me at
% orestisapostolou@yahoo.gr


% ECGanalysis.m must be in the same folder as 106.dat
% it makes a full wavelet analysis with plots for any ECG, in order to
% detect the R peaks

% read the ecgs and their annotations

ann106=rdann('106.dat','atr',1);
[ecg106, Fs, tm] = rdsamp('106.dat',1);


%

ecg = ecg106;
ann = ann106;

N = length(ecg);
x = 1:1:N;
figure
plot(x,ecg106);
hold on
plot(x(ann106(:)),ecg106(ann106(:)),'r*');

%[cA1,cD1] = dwt(ecg106,'coif1');
%A1coif = upcoef('a',cA1,'coif1',1,N);
%D1 = upcoef('d',cD1,'coif1',1,N);

% wavelet decomposition into 2 levels, 'fk4' is the wavelet family

[C,L] = wavedec(ecg,2,'fk4');
[cD1, cD2] = detcoef(C,L,[1,2]);

% D1 is first detail (higher freqs) and D2 is the second

D1 = wrcoef('d',C,L,'fk4',1);
D2 = wrcoef('d',C,L,'fk4',2);

figure
plot(ecg(1:10000));
hold on
plot(D1(1:10000));

% we use D2 because it has lower noise

mx = max(D2);
mn = min(D2);
D2normal(:) = D2(:)-mn;
D2normal(:) = D2normal(:)/(mx-mn);
thres = prctile(D2normal,98);

figure
hist(D2normal,100);

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
ann = nonzeros(ann);
figure
title("o is real, x is predicted");
hold on
plot(x,ecg);
%hold on
%plot(x,D2);
%hold on
%plot(x,y);
hold on
plot(x(locsFinal(:)),ecg(locsFinal(:)),'r*');
hold on
plot(x(ann(:)),ecg(ann(:)),'ko', 'Color', [0 0 0]);

[TP, FP, FN, precision, recall, matchMat] = evaluation(f1,n1,n2,4);

% use a breakpoint to view these metrics

MSE = calculateMSE(ecg, matchMat); 

Rreal = zeros(length(ann)-1,1);
for i = 1:length(Rreal)
    Rreal(i) = ann(i+1) - ann(i);
end

Rpred = zeros(length(locsFinal)-1,1);
for i = 1:length(Rpred)
    Rpred(i) = locsFinal(i+1) - locsFinal(i);
end

figure
plot(Rreal);
hold on
title("Real time interval");
hold off

figure
plot(Rpred);
hold on
title("Predicted time interval");
hold off


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

function MSE = calculateMSE(signal, matchMat)

    [n1,~] = size(matchMat);
    MSE = 0;
    for i = 1: n1
        MSE = MSE + (signal(matchMat(i,1)) - signal(matchMat(i,2))).^2;
    end
    MSE = MSE./n1;


end
