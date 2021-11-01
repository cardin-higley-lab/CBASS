%% Test 1 - Plots the same specrum for different windows in order to see how that scales up
in1WinLen   = [.5 1 2 4 8 16 32];
inNChan     = 15;
hPLT        = nan(size(in1WinLen));
inMaxFreq   = 30;

figure('Position', [50 50 1250 600]); hold on
for iWin = 1:length(in1WinLen)
    in1Idx  = 1:round(sREC.inSampleRate * in1WinLen(iWin));
    cp1FFT  = fft(sREC.db2LFP(inNChan, in1Idx));
    db1Freq = 0:1/in1WinLen(iWin):inMaxFreq;
    hPLT(iWin) = plot(db1Freq, abs(cp1FFT(1:length(db1Freq))));
end
legend(hPLT, cellfun(@num2str, num2cell(in1WinLen), 'UniformOutput', false));
%% Test 2 relation between the zero frequencey and nyquist frequency components
inNItr  = 10000;
db2Edge = nan(inNItr, 2);
in1NFFT = [50 100 200 400 800];
in1MeanMid = nan(size(in1NFFT));

figure('Position', [50 50 1250 600]); hold on
for iNF = 1:length(in1NFFT)
    for iItr = 1:inNItr
        cp1FFT = fft(rand(1, in1NFFT(iNF)));
        db2Edge(iItr, 1) = cp1FFT(1);
        db2Edge(iItr, 2) = cp1FFT((in1NFFT(iNF)/2) + 1);
    end
    % disp(db2Edge);
    in1MeanMid(iNF) = mean(db2Edge(:, 2));
end
plot(in1NFFT, in1MeanMid)
%% Test 3 - relation between the nyquist frequency components when the window is even or uneven
inNItr  = 1000;
db2Edge = nan(inNItr, 3);
inNFFT  = 1000;

for iItr = 1:inNItr
    db1Rand = rand(1,inNFFT + 1);
    cp1FFT  = fft(db1Rand(1:inNFFT - 1));
    db2Edge(iItr, 1) = cp1FFT(inNFFT/2);
    cp1FFT  = fft(db1Rand(1:inNFFT));
    db2Edge(iItr, 2) = cp1FFT((inNFFT/2) + 1);
    cp1FFT  = fft(db1Rand);
    db2Edge(iItr, 3) = cp1FFT((inNFFT/2) + 1);
end

% disp(db2Edge);
figure('Position', [50 50 1250 600])
subplot(1, 3, 1), plot(abs(db2Edge(:, 3)), real(db2Edge(:, 2)), 'x'); title('Abs');
subplot(1, 3, 3), plot(angle(db2Edge(:, 3)), real(db2Edge(:, 2)), 'x'); title('Angle');
subplot(1, 3, 2), plot(real(db2Edge(:, 3)), real(db2Edge(:, 2)), 'x'); title('Real');

figure('Position', [50 50 1250 600])
subplot(1, 2, 1), plot(abs(db2Edge(:, 3)) .* cos(angle(db2Edge(:, 3))), real(db2Edge(:, 2)), 'x'); title('Cos');
subplot(1, 2, 2), plot(abs(db2Edge(:, 3)) .* sin(angle(db2Edge(:, 3))), real(db2Edge(:, 2)), 'x'); title('Sin');

%% Test 4 - fitting procedure on colored noise.
% Generates colored noise
inNSamp     = 1000;
dbClrExp    = 5;
db1Series   = rand(1, inNSamp);
cp1FFT      = fft(db1Series);

% Sets the indices and frequency
in1Idx = 2:inNSamp/2 + 1;
in1Frq = in1Idx - 1; 

% Gets non zero and non Nyquist frequency components and devides them by
% freq
cp1FFT_Sel  = cp1FFT(in1Idx);
% cp1FFT_Clr  = cp1FFT_Sel ./ (sqrt(2 * in1Frq .^ dbClrExp));
cp1FFT_Clr  = cp1FFT_Sel ./ (in1Frq .^ dbClrExp);

% Performs the a less constraining fit
db1Y = 2 .* log(abs(cp1FFT_Clr)) + log(2);
db1X = -log(in1Frq);
dbClrExpHat1 = regress(db1Y', [ones(size(db1Y')) db1X']);

% Performs the a less constraining fit
db1Y = log(abs(cp1FFT_Clr));
db1X = -log(in1Frq);
dbClrExpHat2 = regress(db1Y', [ones(size(db1Y')) db1X']);

% % Plots a figure to see things
figure('Position', [50 50 1250 600]), hold on
hPLT(1) = plot(in1Frq, abs(cp1FFT_Clr));
hPLT(2) = plot(in1Frq, sqrt(exp(dbClrExpHat1(1))./(2 .* in1Frq .^ dbClrExpHat1(2))));
hPLT(3) = plot(in1Frq, exp(dbClrExpHat2(1))./(in1Frq .^ dbClrExpHat2(2)));
set(gca, 'yscale', 'log')
legend(hPLT, {'Data', 'Est1', 'Est2'})
%% test 5 - expected value of white noise
in1NSamp    = [1000 2000 3000 4000 5000];
db1Var      = [1 2 4 8];

db2MeanABS  = nan(length(in1NSamp), length(db1Var));
for iSmp = 1:length(in1NSamp)
    for iMea = 1:length(db1Var)
        db1Series = (rand(1, in1NSamp(iSmp))) .* db1Var(iMea);
        in1Idx = 2:in1NSamp(iSmp)/2 + 1;
        cp1FFT = fft(db1Series);
        db2MeanABS(iSmp, iMea) = mean(abs(cp1FFT(in1Idx)))/sqrt(in1NSamp(iSmp));
    end
end

disp([[nan db1Var]; [in1NSamp' db2MeanABS]])
%% Test 6 - scaling of exponent as a function of window size (i.e. sampling rate)
inNSamp     = 1000;
dbClrExp    = 1;
in1SRate    = [100 200 500 1000 2000 5000];
db1Series   = rand(1, inNSamp);
cp1FFT      = fft(db1Series);
% Sets the indices and frequency
in1Idx      = 2:inNSamp/2 + 1;
in1Frq      = in1Idx - 1; 

% Gets non zero and non Nyquist frequency components and devides them by
% freq
cp1FFT_Sel  = cp1FFT(in1Idx);
cp1FFT_Clr  = cp1FFT_Sel ./ (in1Frq .^ dbClrExp);
db1Y        = log(abs(cp1FFT_Clr));

db2Beta = nan(2, length(in1SRate));
for iRate = 1:length(in1SRate)
    in1Frq_Fit          = in1Frq .* in1SRate(iRate) ./ inNSamp;
    db1X                = -log(in1Frq_Fit);
    db2Beta(:, iRate)   = regress(db1Y', [ones(size(db1Y')) db1X']);
end

%Plots some figures
hFIG = figure('Position', [50 50 1250 600]);
subplot(1, 3, 1);
clear hPLT
plot(in1Frq, abs(cp1FFT_Clr)); hold on
for iPlt = 1:length(in1SRate)
    in1Frq_Fit  = in1Frq .* in1SRate(iPlt) ./ inNSamp;
    hPLT(iPlt)  = plot(in1Frq + iPlt * 0.1, exp(db2Beta(1, iPlt))./(in1Frq_Fit .^ db2Beta(2, iPlt)));
end
set(gca, 'yscale', 'log')
% legend(hPLT, cellfun(@num2str, num2cell(in1SRate), 'UniformOutput', false));
legend(cellfun(@num2str, num2cell(in1SRate), 'UniformOutput', false));
title('Fits')
subplot(1, 3, 2)
plot(in1SRate, db2Beta(1, :), 'o--'); title('Intercept')
subplot(1, 3, 3)
plot(in1SRate, db2Beta(2, :), 'o--'); title('Exponent')